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
        tau,
        tau_limit,
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

    void EmissionParticles(const amrex::MultiFab& matter, const amrex::Real n_nu_packet, const amrex::Real nu_Energy_MeV,  const amrex::Real dtdE3_3dOmegadx3, const amrex::Real curr_time_s, const int simtype);

    void AbsorptionParticles();

};

#endif