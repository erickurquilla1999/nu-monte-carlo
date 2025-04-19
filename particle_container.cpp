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
MCParticleContainer::InitParticlesTest1(const int num_par_test_1)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    amrex::Print() << "Cell size (dx): " << dx[0] << ", " << dx[1] << ", " << dx[2] << "\n";
    const Real* plo = geom.ProbLo();
    amrex::Print() << "ProbLo (plo): " << plo[0] << ", " << plo[1] << ", " << plo[2] << "\n";

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};

        const int grid_id = mfi.index();
        amrex::Print() << "Grid ID: " << grid_id << "\n";
        const int tile_id = mfi.LocalTileIndex();
        amrex::Print() << "Tile ID: " << tile_id << "\n";
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            amrex::Print() << "IntVect iv: (" << iv[0] << ", " << iv[1] << ", " << iv[2] << ")\n";

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
                    p.rdata(RealData::time_s) = 0.0;
                    p.rdata(RealData::N) = 1.0;

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
    const int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        amrex::Print() << "Processing particle tile: " << pti.index() << ", " << pti.LocalTileIndex() << "\n";
        const int np  = pti.numParticles();
        amrex::Print() << "Number of particles in tile: " << np << "\n";
        ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);

        amrex::ParallelFor (np, [=] AMREX_GPU_DEVICE (int i) {

            ParticleType& p = pstruct[i];

            printf("Particle position pos: (%f, %f, %f)\n",
                   p.pos(0),
                   p.pos(1),
                   p.pos(2));
            printf("Particle position rea: (%f, %f, %f)\n",
                    p.rdata(RealData::x),
                    p.rdata(RealData::y),
                    p.rdata(RealData::z));

            printf("p.id(): %d\n", p.id());
            printf("p.cpu(): %d\n", p.cpu());
        });
    }
}

void
MCParticleContainer::MoveParticles(const amrex::Real dt)
{
    const int lev = 0;
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        amrex::Print() << "Processing particle tile: " << pti.index() << ", " << pti.LocalTileIndex() << "\n";
        const int np  = pti.numParticles();
        amrex::Print() << "Number of particles in tile: " << np << "\n";
        ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);

        amrex::ParallelFor (np, [=] AMREX_GPU_DEVICE (int i) {

            ParticleType& p = pstruct[i];

            amrex::Real c_cm_s = 3.0e10; // speed of light in cm/s

            p.rdata(RealData::x) += p.rdata(RealData::phatx) * dt * c_cm_s;
            p.rdata(RealData::y) += p.rdata(RealData::phaty) * dt * c_cm_s;
            p.rdata(RealData::z) += p.rdata(RealData::phatz) * dt * c_cm_s;
            p.rdata(RealData::time_s) += dt;


        });
    }
}