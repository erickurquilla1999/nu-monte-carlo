#include <AMReX_MultiFab.H>

#include "matter.h"

void 
init_matter(amrex::MultiFab& matter)
{
    for (amrex::MFIter mfi(matter); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::Array4<amrex::Real>& mf_array = matter.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            if (i == 0 && j == 0 && k == 0) {
                mf_array(i,j,k,MatterData::IMFP_cm) = 0.0; // 1/cm
            } else {
                mf_array(i,j,k,MatterData::IMFP_cm) = 1.0/5.5; // 1/cm
            }
            #ifdef DEBUG
                printf("mf_array(%d,%d,%d,%d) = %f\n", i, j, k, MatterData::rho_g_ccm, mf_array(i,j,k,MatterData::rho_g_ccm)); // Print the value
            #endif
        });
    }
}
