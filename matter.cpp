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

void
compute_nu_per_MC_particles(amrex::MultiFab& matter, int n_mc_particles, int& n_nu_per_mc_particles, const amrex::Real nu_Energy_MeV)
{

    amrex::ReduceOps< amrex::ReduceOpSum, amrex::ReduceOpMin, amrex::ReduceOpMax > reduce_op;
    amrex::ReduceData< amrex::Real , amrex::Real, amrex::Real > reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (amrex::MFIter mfi(matter); mfi.isValid(); ++mfi) {

        const amrex::Box& bx = mfi.validbox();
        auto const& matter_multifab = matter.array(mfi);

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {

            amrex::Real feq = 1.0 / ( 1.0 + exp( ( nu_Energy_MeV - matter_multifab(i,j,k,MatterData::chemical_potential_MeV) ) / matter_multifab(i,j,k,MatterData::T_MeV) ) );

            return {feq, matter_multifab(i,j,k,MatterData::IMFP_cm), matter_multifab(i,j,k,MatterData::IMFP_cm)};
        });
    }

    // extract the reduced values from the combined reduced data structure
    auto rv = reduce_data.value();
    amrex::Real sumfeq = amrex::get<0>(rv);
    amrex::Real minIMFPcm   = amrex::get<1>(rv);
    amrex::Real maxIMFPcm   = amrex::get<2>(rv);

    // reduce across MPI ranks
    amrex::ParallelDescriptor::ReduceRealMin(sumfeq);
    amrex::ParallelDescriptor::ReduceRealMin(minIMFPcm);
    amrex::ParallelDescriptor::ReduceRealMin(maxIMFPcm);

}
