#ifndef MATTER_H
#define MATTER_H

struct MatterData
{
    enum {
        rho_g_ccm,
        ye,
        T_MeV,
        IMFP_cm,
        chemical_potential_MeV,
        ncomps
    };
};

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

void init_matter(amrex::MultiFab& matter);

void compute_nu_per_MC_particles(amrex::MultiFab& matter, const int n_mc_particles, amrex::Real& n_nu_per_mc_particles, const amrex::Real nu_Energy_MeV, const amrex::Real dtdE3_3dOmegadx3);

#endif