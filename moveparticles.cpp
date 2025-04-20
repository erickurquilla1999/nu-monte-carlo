
#include "AMReX_MultiFab.H"
#include "AMReX_ParticleContainer.H"
#include <AMReX.H>

# include "moveparticles.h"
#include "particle_container.h"
#include "matter.h"

void MoveParticlesMC(MCParticleContainer& particles, const amrex::MultiFab& state, const amrex::Geometry& geom, const amrex::Real dt)
{

    const amrex::Real* dx = geom.CellSize();
    const amrex::Real* plo = geom.ProbLo();
    amrex::Real dx_local[3] = {dx[0], dx[1], dx[2]};
    amrex::Real plo_local[3] = {plo[0], plo[1], plo[2]};

    amrex::MeshToParticle(particles, state, 0,
    [=] AMREX_GPU_DEVICE (MCParticleContainer::ParticleType& p,
                          amrex::Array4<const amrex::Real> const& sarr)
    {
        amrex::Real imfp_cm = sarr(p.idata(IntData::i), p.idata(IntData::j), p.idata(IntData::k), MatterData::IMFP_cm);

        amrex::Real c_cm_s = 3.0e10; // speed of light in cm/s

        p.pos(0) += p.rdata(RealData::phatx) * dt * c_cm_s;
        p.pos(0) += p.rdata(RealData::phaty) * dt * c_cm_s;
        p.pos(0) += p.rdata(RealData::phatz) * dt * c_cm_s;

        p.rdata(RealData::x) += p.rdata(RealData::phatx) * dt * c_cm_s;
        p.rdata(RealData::y) += p.rdata(RealData::phaty) * dt * c_cm_s;
        p.rdata(RealData::z) += p.rdata(RealData::phatz) * dt * c_cm_s;
        p.rdata(RealData::time_s) += dt;

        p.rdata(RealData::tau) += c_cm_s * dt * imfp_cm;

        if (p.rdata(RealData::tau) > p.rdata(RealData::tau_limit)) {
            p.pos(0) = -1.0;
            p.pos(1) = -1.0;
            p.pos(2) = -1.0;
        }
        
    });
}
