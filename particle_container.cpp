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
MCParticleContainer::InitParticlesTest1(const int nbox, const int num_par_test_1)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    amrex::Print() << "Cell size (dx): " << dx[0] << ", " << dx[1] << ", " << dx[2] << "\n";
    const Real* plo = geom.ProbLo();
    amrex::Print() << "ProbLo (plo): " << plo[0] << ", " << plo[1] << ", " << plo[2] << "\n";

    const int num_par_per_tile = num_par_test_1 / nbox;
    amrex::Print() << "Number of particles per tile: " << num_par_per_tile << "\n";

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};

        const int grid_id = mfi.index();
        amrex::Print() << "Grid ID: " << grid_id << "\n";
        const int tile_id = mfi.LocalTileIndex();
        amrex::Print() << "Tile ID: " << tile_id << "\n";

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        IntVect init_indx_cell_in_tile_box = tile_box.smallEnd();

        amrex::Print() << "IntVect init_indx_cell_in_tile_box: (" << init_indx_cell_in_tile_box[0] << ", " << init_indx_cell_in_tile_box[1] << ", " << init_indx_cell_in_tile_box[2] << ")\n";

        for (int i_part=0; i_part<num_par_per_tile;i_part++) {

            ParticleType p;

            p.id()  = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            p.pos(0) = plo[0] + (init_indx_cell_in_tile_box[0] + 0.5) * dx[0]; 
            p.pos(1) = plo[1] + (init_indx_cell_in_tile_box[1] + 0.5) * dx[1];
            p.pos(2) = plo[2] + (init_indx_cell_in_tile_box[2] + 0.5) * dx[2];

            p.rdata(RealData::x) = plo[0] + dx[0] * 0.5;
            p.rdata(RealData::y) = plo[1] + dx[1] * 0.5;
            p.rdata(RealData::z) = plo[2] + dx[2] * 0.5;
            p.rdata(RealData::E_MeV) = 1.0;
            p.rdata(RealData::phatx) = 1.0;
            p.rdata(RealData::phaty) = 0.0;
            p.rdata(RealData::phatz) = 0.0;
            p.rdata(RealData::time) = 0.0;
            p.rdata(RealData::N) = 1.0;

            // AMREX_ASSERT(this->Index(p, lev) == iv);

            particle_tile.push_back(p);

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