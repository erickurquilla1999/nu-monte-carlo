#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "parameters.h"

void evolve()
{
    amrex::Print() << "AMReX version " << amrex::Version() << "\n";

    // Get the simulation parameters
    Parameters_struct params;
    get_parameters(params);

    amrex::Print() << "Input parameters:\n";
    amrex::Print() << "  input_n_cell_x: " << params.input_n_cell_x << "\n";
    amrex::Print() << "  input_n_cell_y: " << params.input_n_cell_y << "\n";
    amrex::Print() << "  input_n_cell_z: " << params.input_n_cell_z << "\n";
    amrex::Print() << "  input_max_grid_size_x: " << params.input_max_grid_size_x << "\n";
    amrex::Print() << "  input_max_grid_size_y: " << params.input_max_grid_size_y << "\n";
    amrex::Print() << "  input_max_grid_size_z: " << params.input_max_grid_size_z << "\n";
    amrex::Print() << "  input_physical_domain_size_x_cm: " << params.input_physical_domain_size_x_cm << "\n";
    amrex::Print() << "  input_physical_domain_size_y_cm: " << params.input_physical_domain_size_y_cm << "\n";
    amrex::Print() << "  input_physical_domain_size_z_cm: " << params.input_physical_domain_size_z_cm << "\n";

    // Set up the simulation domain
    amrex::IntVect domain_low_cells(0,0,0);
    amrex::IntVect domain_high_cells(params.input_n_cell_x, params.input_n_cell_y, params.input_n_cell_z);

    // Create the domain box
    amrex::Box domain_box(domain_low_cells, domain_high_cells);
    amrex::BoxArray domain_box_array(domain_box);

    // Set the number of subdivisions for each box
    amrex::IntVect max_grid_size(params.input_max_grid_size_x, params.input_max_grid_size_y, params.input_max_grid_size_z);
    domain_box_array.maxSize(max_grid_size); // Set the maximum size of each box

    // Create the distribution mapping
    amrex::DistributionMapping distribution_mapping(domain_box_array);

    int number_componets = 1;
    int number_ghost_cells = 0;

    amrex::MultiFab mf_Ye(domain_box_array, distribution_mapping, number_componets, number_ghost_cells);
    amrex::MultiFab mf_rho_g_ccm(domain_box_array, distribution_mapping, number_componets, number_ghost_cells);
    amrex::MultiFab mf_T_MeV(domain_box_array, distribution_mapping, number_componets, number_ghost_cells);
    amrex::MultiFab mf_absorption_imfp_cm(domain_box_array, distribution_mapping, number_componets, number_ghost_cells);

    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(params.input_physical_domain_size_x_cm, params.input_physical_domain_size_y_cm, params.input_physical_domain_size_z_cm)});

    amrex::Geometry geom(domain_box, &real_box);
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> cell_size_cm = geom.CellSizeArray();

    for (amrex::MFIter mfi(mf_absorption_imfp_cm); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const amrex::Array4<amrex::Real>& mf_array_absorption_imfp_cm = mf_absorption_imfp_cm.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            amrex::Real x = (i + 0.5) * cell_size_cm[0];
            amrex::Real y = (j + 0.5) * cell_size_cm[1];
            amrex::Real z = (k + 0.5) * cell_size_cm[2];

            mf_array_absorption_imfp_cm(i,j,k) = x; // Set the value of the MultiFab to 1.0
        });
    }

    WriteSingleLevelPlotfile("plt001", mf_absorption_imfp_cm, {"absorption_imfp_cm"}, geom, 0.0, 0);
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    evolve();
    amrex::Finalize();
}