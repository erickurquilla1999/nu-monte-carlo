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
                p.pos(0)             = plo[0] + (i + 0.5*(rand+1.0))*dx[0];
                p.rdata(RealData::x) = plo[0] + (i + 0.5*(rand+1.0))*dx[0];

                symmetric_uniform(&rand, engine);
                p.pos(1)             = plo[1] + (j + 0.5*(rand+1.0))*dx[1];
                p.rdata(RealData::y) = plo[1] + (j + 0.5*(rand+1.0))*dx[1];

                symmetric_uniform(&rand, engine);
                p.pos(2)             = plo[2] + (k + 0.5*(rand+1.0))*dx[2];
                p.rdata(RealData::z) = plo[2] + (k + 0.5*(rand+1.0))*dx[2];

                // Set particle momentum

                symmetric_uniform(&rand, engine);
                amrex::Real phi = ((rand+1.0)*0.5)*2.0*MathConst::pi;
                symmetric_uniform(&rand, engine);
                amrex::Real cos_theta = rand;
                amrex::Real sin_theta = sqrt(1 - cos_theta*cos_theta);

                p.rdata(RealData::phatx) = sin_theta * cos(phi);
                p.rdata(RealData::phaty) = sin_theta * sin(phi);
                p.rdata(RealData::phatz) = cos_theta;

                // Test values
                // p.rdata(RealData::phatx) = 1.0;
                // p.rdata(RealData::phaty) = 0.0;
                // p.rdata(RealData::phatz) = 0.0;

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