#ifndef MOVE_PARTICLES_H
#define MOVE_PARTICLES_H

#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>

#include "particle_container.h"

struct NeutrinoData
{
    enum {
        n_invcm3,
        fx_invcm3,
        fy_invcm3,
        fz_invcm3,
        ncomps
    };
};

void MoveParticlesMC(MCParticleContainer& particles, const amrex::MultiFab& state, const amrex::Geometry& geom, const amrex::Real dt);

void compute_nu_n_and_f(MCParticleContainer& particles, const amrex::Geometry& geom, amrex::MultiFab& nu_n_and_f);

#endif