#ifndef MATTER_H
#define MATTER_H

#include <AMReX_MultiFab.H>

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

void init_matter(amrex::MultiFab& matter);

#endif