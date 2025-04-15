#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H

#include <AMReX_Particles.H>

struct RealData
{
    enum {
        time = 0,
        x,
        y,
        z,
        pupx,
        pupy,
        pupz,
        pupt,
        ncomps
    };
};

struct IntData
{
    enum {
        i,
        j,
        k,
        ncomps
    };
};

class MCParticleContainer
    : public amrex::ParticleContainer<RealData::ncomps, IntData::ncomps>
{
public:

    using MyParIter = amrex::ParIter<RealData::ncomps, IntData::ncomps>;

    MCParticleContainer(const amrex::Geometry            & a_geom,
                        const amrex::DistributionMapping & a_dmap,
                        const amrex::BoxArray            & a_ba);

    void InitParticles();

    void LoopParticlesPrint();

};

#endif