#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H

#include <AMReX_Particles.H>

struct RealData
{
    enum {
        time_s = 0,
        x,
        y,
        z,
        phatx,
        phaty,
        phatz,
        E_MeV,
        N,
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

    void InsertParticles(const int num_par_test_1);

    void LoopParticlesPrint();

    void MoveParticles(const amrex::Real dt);

};

#endif