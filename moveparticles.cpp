
#include "AMReX_MultiFab.H"
#include "AMReX_ParticleContainer.H"
#include <AMReX.H>

# include "moveparticles.h"
#include "particle_container.h"

void MoveParticlesMC(MCParticleContainer& particles, const amrex::MultiFab& state, const amrex::Geometry& geom, const amrex::Real dt)
{

    const amrex::Real* dx = geom.CellSize();
    const amrex::Real* plo = geom.ProbLo();
    amrex::Real dx_local[3] = {dx[0], dx[1], dx[2]};
    amrex::Real plo_local[3] = {plo[0], plo[1], plo[2]};

    amrex::Print() << "dx: " << dx_local[0] << ", " << dx_local[1] << ", " << dx_local[2] << "\n";
    amrex::Print() << "plo: " << plo_local[0] << ", " << plo_local[1] << ", " << plo_local[2] << "\n";
    
    amrex::MeshToParticle(particles, state, 0,
    [=] AMREX_GPU_DEVICE (MCParticleContainer::ParticleType& p,
                          amrex::Array4<const amrex::Real> const& sarr)
    {

        
        
        // // The following variables contains temperature, electron fraction, and density interpolated from grid quantities to particle positions
        // Real T_pp = 0; // erg
        // Real Ye_pp = 0;
        // Real rho_pp = 0; // g/ccm

        // for (int k = sz.first(); k <= sz.last(); ++k) {
        //     for (int j = sy.first(); j <= sy.last(); ++j) {
        //         for (int i = sx.first(); i <= sx.last(); ++i) {
        //             #include "generated_files/Evolve.cpp_interpolate_from_mesh_fill"
        //         }
        //     }
        // }
        
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
