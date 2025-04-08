#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include "parameters.h"

void evolve()
{
        amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";

        Parameters_struct params;
        get_parameters(params);

        amrex::IntVect domain_low_cells(0,0,0);
        amrex::IntVect domain_high_cells(params.input_n_cell_x, params.input_n_cell_y, params.input_n_cell_z);

        amrex::Print() << "Input parameters:\n";
        amrex::Print() << "  input_n_cell_x: " << params.input_n_cell_x << "\n";

        amrex::Box domain_box(domain_low_cells, domain_high_cells);
        amrex::BoxArray domain_box_array(domain_box);

        amrex::IntVect max_grid_size(params.input_max_grid_size_x, params.input_max_grid_size_y, params.input_max_grid_size_z); // 
        domain_box_array.maxSize(max_grid_size); // Set the maximum size of each box

        amrex::DistributionMapping distribution_mapping(domain_box_array);

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    evolve();
    amrex::Finalize();
}