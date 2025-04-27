#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "particle_container.h"

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