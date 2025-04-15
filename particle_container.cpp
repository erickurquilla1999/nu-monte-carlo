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
MCParticleContainer::InitParticles()
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real* dx = geom.CellSize();
    const Real* plo = geom.ProbLo();

    const int num_ppc = 10;

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
        {
            for (int i_part=0; i_part<num_ppc;i_part++) {
                
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();

                p.pos(0) = 1.0;
                p.pos(1) = 1.0;
                p.pos(2) = 1.0;

                // p.rdata(RealData::vx) = u[0];
                // p.rdata(RealData::vy) = u[1];
                // p.rdata(RealData::vz) = u[2];

                // AMREX_ASSERT(this->Index(p, lev) == iv);

                particle_tile.push_back(p);
            }
        }
    }


}