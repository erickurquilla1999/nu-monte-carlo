
#include "AMReX_MultiFab.H"
#include "AMReX_ParticleContainer.H"
#include <AMReX.H>

# include "particle_interactions.h"
#include "particle_container.h"
#include "matter.h"
#include "constant.h"

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
        int i_indx = amrex::Math::floor((p.pos(0) - plo_local[0]) / dx_local[0]);
        int j_indx = amrex::Math::floor((p.pos(1) - plo_local[1]) / dx_local[1]);
        int k_indx = amrex::Math::floor((p.pos(2) - plo_local[2]) / dx_local[2]);

        amrex::Real imfp_cm = sarr(i_indx, j_indx, k_indx, MatterData::IMFP_cm);

        p.pos(0) += p.rdata(RealData::phatx) * dt * PhysConst::c;
        p.pos(1) += p.rdata(RealData::phaty) * dt * PhysConst::c;
        p.pos(2) += p.rdata(RealData::phatz) * dt * PhysConst::c;

        p.rdata(RealData::x) += p.rdata(RealData::phatx) * dt * PhysConst::c;
        p.rdata(RealData::y) += p.rdata(RealData::phaty) * dt * PhysConst::c;
        p.rdata(RealData::z) += p.rdata(RealData::phatz) * dt * PhysConst::c;
        p.rdata(RealData::time_s) += dt;

        p.rdata(RealData::tau) += PhysConst::c * dt * imfp_cm;

        if (p.rdata(RealData::tau) > p.rdata(RealData::tau_limit)) {
            p.pos(0) = -1.0;
            p.pos(1) = -1.0;
            p.pos(2) = -1.0;
        }
        
    });
}

void compute_nu_n_and_f(MCParticleContainer& particles, const amrex::Geometry& geom, amrex::MultiFab& nu_n_and_f)
{
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const amrex::Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];

    nu_n_and_f.setVal(0.0);

    amrex::ParticleToMesh(particles, nu_n_and_f, 0,
    [=] AMREX_GPU_DEVICE (const MCParticleContainer::ParticleType& p,
                          amrex::Array4<amrex::Real> const& sarr)
    {
        int i_indx = amrex::Math::floor((p.pos(0) - plo[0]) * dxi[0]);
        int j_indx = amrex::Math::floor((p.pos(1) - plo[1]) * dxi[1]);
        int k_indx = amrex::Math::floor((p.pos(2) - plo[2]) * dxi[2]);

        amrex::Gpu::Atomic::AddNoRet(&sarr(i_indx,j_indx,k_indx,NeutrinoData::n_invcm3), p.rdata(RealData::N)*inv_cell_volume);
        amrex::Gpu::Atomic::AddNoRet(&sarr(i_indx,j_indx,k_indx,NeutrinoData::fx_invcm3), p.rdata(RealData::N)*p.rdata(RealData::phatx)*inv_cell_volume);
        amrex::Gpu::Atomic::AddNoRet(&sarr(i_indx,j_indx,k_indx,NeutrinoData::fy_invcm3), p.rdata(RealData::N)*p.rdata(RealData::phaty)*inv_cell_volume);
        amrex::Gpu::Atomic::AddNoRet(&sarr(i_indx,j_indx,k_indx,NeutrinoData::fz_invcm3), p.rdata(RealData::N)*p.rdata(RealData::phatz)*inv_cell_volume);
    });
}