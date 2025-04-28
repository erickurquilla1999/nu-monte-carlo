#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "particle_container.h"
#include "constant.h"
#include "matter.h"

using namespace amrex;

MCParticleContainer::
MCParticleContainer(const Geometry            & a_geom,
                    const DistributionMapping & a_dmap,
                    const BoxArray            & a_ba)
    : ParticleContainer<RealData::ncomps, IntData::ncomps> (a_geom, a_dmap, a_ba)
{
}

void
MCParticleContainer::InsertParticles(const int num_par_test_1, amrex::Real current_time_s)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();

    #ifdef DEBUG
        amrex::Print() << "Cell size (dx): " << dx[0] << ", " << dx[1] << ", " << dx[2] << "\n";
        amrex::Print() << "ProbLo (plo): " << plo[0] << ", " << plo[1] << ", " << plo[2] << "\n";
    #endif


    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();

        #ifdef DEBUG
            amrex::Print() << "Grid ID: " << grid_id << "\n";
            amrex::Print() << "Tile ID: " << tile_id << "\n";
        #endif

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            #ifdef DEBUG
                amrex::Print() << "IntVect iv: (" << iv[0] << ", " << iv[1] << ", " << iv[2] << ")\n";
            #endif

            if (iv[0] == 0 && iv[1] == 0 && iv[2] == 0) {

                for (int i_part=0; i_part<num_par_test_1;i_part++) {

                    ParticleType p;

                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();

                    p.pos(0) = plo[0] + dx[0] * 0.5;
                    p.pos(1) = plo[1] + dx[1] * 0.5;
                    p.pos(2) = plo[2] + dx[2] * 0.5;

                    p.rdata(RealData::x) = plo[0] + dx[0] * 0.5;
                    p.rdata(RealData::y) = plo[1] + dx[1] * 0.5;
                    p.rdata(RealData::z) = plo[2] + dx[2] * 0.5;
                    p.rdata(RealData::E_MeV) = 1.0;
                    p.rdata(RealData::phatx) = 1.0;
                    p.rdata(RealData::phaty) = 0.0;
                    p.rdata(RealData::phatz) = 0.0;
                    p.rdata(RealData::time_s) = current_time_s;
                    p.rdata(RealData::N) = 1.0;
                    p.rdata(RealData::tau) = 0.0;
                    p.rdata(RealData::tau_limit) = -std::log(amrex::Random());
                    
                    p.idata(IntData::i) = iv[0];
                    p.idata(IntData::j) = iv[1];
                    p.idata(IntData::k) = iv[2];

                    // AMREX_ASSERT(this->Index(p, lev) == iv);

                    particle_tile.push_back(p);

                }
            }
        }
    }
}

void
MCParticleContainer::LoopParticlesPrint()
{
 
    printf("time_s x y z phatx phaty phatz E_MeV N tau tau_limit i j k\n");

    const int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        // amrex::Print() << "Processing particle tile: " << pti.index() << ", " << pti.LocalTileIndex() << "\n";
        const int np  = pti.numParticles();
        // amrex::Print() << "Number of particles in tile: " << np << "\n";
        ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);

        amrex::ParallelFor (np, [=] AMREX_GPU_DEVICE (int i) {

            ParticleType& p = pstruct[i];

            // printf("Particle position pos: (%f, %f, %f)\n",
            //        p.pos(0),
            //        p.pos(1),
            //        p.pos(2));
            // printf("p.id(): %d\n", static_cast<int>(p.id()));
            // printf("p.cpu(): %d\n", static_cast<int>(p.cpu()));

            printf("%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %d %d %d\n",
                p.rdata(RealData::time_s),
                p.rdata(RealData::x),
                p.rdata(RealData::y),
                p.rdata(RealData::z),
                p.rdata(RealData::phatx),
                p.rdata(RealData::phaty),
                p.rdata(RealData::phatz),
                p.rdata(RealData::E_MeV),
                p.rdata(RealData::N),
                p.rdata(RealData::tau),
                p.rdata(RealData::tau_limit),
                static_cast<int>(p.idata(IntData::i)),
                static_cast<int>(p.idata(IntData::j)),
                static_cast<int>(p.idata(IntData::k))
                );


        });
    }
}

void
MCParticleContainer::MoveParticles(const amrex::Real dt)
{
    const int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {

        const int np  = pti.numParticles();
        ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);

        amrex::ParallelFor (np, [=] AMREX_GPU_DEVICE (int i) {

            ParticleType& p = pstruct[i];

            amrex::Real c_cm_s = 3.0e10; // speed of light in cm/s

            p.pos(0) += p.rdata(RealData::phatx) * dt * c_cm_s;
            p.pos(0) += p.rdata(RealData::phaty) * dt * c_cm_s;
            p.pos(0) += p.rdata(RealData::phatz) * dt * c_cm_s;

            p.rdata(RealData::x) += p.rdata(RealData::phatx) * dt * c_cm_s;
            p.rdata(RealData::y) += p.rdata(RealData::phaty) * dt * c_cm_s;
            p.rdata(RealData::z) += p.rdata(RealData::phatz) * dt * c_cm_s;
            p.rdata(RealData::time_s) += dt;


        });
    }
}

void
MCParticleContainer::UpdateCellIndex()
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();
    Real dx_local[3] = {dx[0], dx[1], dx[2]};
    Real plo_local[3] = {plo[0], plo[1], plo[2]};

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {

        const int np  = pti.numParticles();
        ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);

        amrex::ParallelFor (np, [=] AMREX_GPU_DEVICE (int i) {
            ParticleType& p = pstruct[i];

            p.idata(IntData::i) = amrex::Math::floor((p.pos(0) - plo_local[0]) / dx_local[0]);
            p.idata(IntData::j) = static_cast<int>( (p.pos(1) - plo_local[1]) / dx_local[1] );
            p.idata(IntData::k) = static_cast<int>( (p.pos(2) - plo_local[2]) / dx_local[2] );

        });
    }
}





namespace
{
    AMREX_GPU_HOST_DEVICE void symmetric_uniform(Real* Usymmetric, amrex::RandomEngine const& engine)
    {
        *Usymmetric = 2. * (amrex::Random(engine)-0.5);
    }
}





void
MCParticleContainer::
EmissionParticles(const amrex::MultiFab& matter, const amrex::Real n_nu_packet, const amrex::Real nu_Energy_MeV,  const amrex::Real dtdE3_3dOmegadx3, const amrex::Real curr_time_s)
{

    const int lev = 0;
    const auto dx = Geom(lev).CellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();
    const auto& a_bounds = Geom(lev).ProbDomain();


    // determine the number of directions per location

    const Real scale_fac = dx[0]*dx[1]*dx[2];
    printf("scale_fac: %f\n", scale_fac);
    printf("dx: %f %f %f\n", dx[0], dx[1], dx[2]);
    printf("plo: %f %f %f\n", plo[0], plo[1], plo[2]);


    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();

        const auto lo = amrex::lbound(tile_box);
        const auto hi = amrex::ubound(tile_box);

        amrex::Print() << "hi.x: " << hi.x << ", hi.y: " << hi.y << ", hi.z: " << hi.z << "\n";
        amrex::Print() << "lo.x: " << lo.x << ", lo.y: " << lo.y << ", lo.z: " << lo.z << "\n";

        Gpu::ManagedVector<unsigned int> counts(tile_box.numPts(), 0);
        unsigned int* pcount = counts.dataPtr();

        Gpu::ManagedVector<unsigned int> offsets(tile_box.numPts());
        unsigned int* poffset = offsets.dataPtr();

        auto const& matter_multifab = matter.array(mfi);

        // Determine how many particles to add to the particle tile per cell
        amrex::ParallelFor(tile_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            amrex::Real feq = 1.0 / ( 1.0 + exp( ( nu_Energy_MeV - matter_multifab(i,j,k,MatterData::chemical_potential_MeV) ) / matter_multifab(i,j,k,MatterData::T_MeV) ) );
            amrex::Real deltaN = ( 1.0 / ( PhysConst::c2 * PhysConst::hbar * PhysConst::hbar * PhysConst::hbar ) ) * dtdE3_3dOmegadx3 * matter_multifab(i,j,k,MatterData::IMFP_cm) * feq;
            int num_to_add = static_cast<int>( deltaN / n_nu_packet );

            int ix = i - lo.x;
            int iy = j - lo.y;
            int iz = k - lo.z;
            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;
            unsigned int uix = amrex::min(nx-1,amrex::max(0,ix));
            unsigned int uiy = amrex::min(ny-1,amrex::max(0,iy));
            unsigned int uiz = amrex::min(nz-1,amrex::max(0,iz));
            unsigned int cellid = (uix * ny + uiy) * nz + uiz;
            pcount[cellid] += num_to_add;
        });

        // Determine total number of particles to add to the particle tile
        Gpu::inclusive_scan(counts.begin(), counts.end(), offsets.begin());

        int num_to_add = offsets[tile_box.numPts()-1];
        amrex::Print() << "Number of particles to add: " << num_to_add << "\n";
        if (num_to_add == 0) continue;

        // this will be the particle ID for the first new particle in the tile
        long new_pid;
        ParticleType* pstruct;

        // #ifdef _OPENMP
        // #pragma omp critical
        // #endif
        // {
            auto& particles = GetParticles(lev);
            auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

            // Resize the particle container
            auto old_size = particle_tile.GetArrayOfStructs().size();
            printf("old_size: %lu\n", old_size);
            auto new_size = old_size + num_to_add;
            printf("new_size: %lu\n", new_size);

            particle_tile.resize(new_size);

            // get the next particle ID
            new_pid = ParticleType::NextID();

            // set the starting particle ID for the next tile of particles
            ParticleType::NextID(new_pid + num_to_add);

            pstruct = particle_tile.GetArrayOfStructs()().data();
        // }

        int procID = ParallelDescriptor::MyProc();

        //===============================================//
        // Initialize particle data in the particle tile //
        //===============================================//
        amrex::ParallelForRNG(tile_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
            int ix = i - lo.x;
            int iy = j - lo.y;
            int iz = k - lo.z;
            int nx = hi.x-lo.x+1;
            int ny = hi.y-lo.y+1;
            int nz = hi.z-lo.z+1;
            unsigned int uix = amrex::min(nx-1,amrex::max(0,ix));
            unsigned int uiy = amrex::min(ny-1,amrex::max(0,iy));
            unsigned int uiz = amrex::min(nz-1,amrex::max(0,iz));
            unsigned int cellid = (uix * ny + uiy) * nz + uiz;

            printf("pcount[%u]: %u\n", cellid, pcount[cellid]);

            for (int newparthiscell=0; newparthiscell<pcount[cellid];newparthiscell++)
            {

                // Get the Particle data corresponding to our particle index in pidx
                const int pidx = poffset[cellid] - poffset[0] + old_size + newparthiscell;
                printf("pidx: %d\n", pidx);

                ParticleType& p = pstruct[pidx];

                // Set particle ID using the ID for the first of the new particles in this tile
                // plus our zero-based particle index
                p.id()   = new_pid + pidx;
                printf("p.id(): %d\n", static_cast<int>(p.id()));

                // Set CPU ID
                p.cpu()  = procID;

                // set particle time
                p.rdata(RealData::time_s) = curr_time_s;

                // Set particle position
                Real rand;

                symmetric_uniform(&rand, engine);
                p.pos(0)             = plo[0] + (i + 0.5*rand)*dx[0];
                p.rdata(RealData::x) = plo[0] + (i + 0.5*rand)*dx[0];

                symmetric_uniform(&rand, engine);
                p.pos(1)             = plo[1] + (j + 0.5*rand)*dx[1];
                p.rdata(RealData::y) = plo[1] + (j + 0.5*rand)*dx[1];

                symmetric_uniform(&rand, engine);
                p.pos(2)             = plo[2] + (k + 0.5*rand)*dx[2];
                p.rdata(RealData::z) = plo[2] + (k + 0.5*rand)*dx[2];

                // Set particle momentum

                symmetric_uniform(&rand, engine);
                amrex::Real phi = ((rand+1.0)*0.5)*2.0*MathConst::pi;
                symmetric_uniform(&rand, engine);
                amrex::Real cos_theta = rand;
                amrex::Real sin_theta = sqrt(1 - cos_theta*cos_theta);

                p.rdata(RealData::phatx) = sin_theta * cos(phi);
                p.rdata(RealData::phaty) = sin_theta * sin(phi);
                p.rdata(RealData::phatz) = cos_theta;

                // Set particle energy
                p.rdata(RealData::E_MeV) = nu_Energy_MeV;

                // Set number of physical particles in MC particle
                p.rdata(RealData::N) = n_nu_packet;

                // Set particle optical depth
                p.rdata(RealData::tau) = 0.0;
                symmetric_uniform(&rand, engine);
                p.rdata(RealData::tau_limit) = -std::log(((rand+1.0)*0.5));
            }

        }); // loop over grid cells


    } // loop over multifabs

} // InitParticles()