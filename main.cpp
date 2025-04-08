#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

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

        // Set up the simulation domain
        amrex::IntVect domain_low_cells(0,0,0);
        amrex::IntVect domain_high_cells(params.input_n_cell_x, params.input_n_cell_y, params.input_n_cell_z);

        // Create the domain box
        amrex::Box domain_box(domain_low_cells, domain_high_cells);
        amrex::BoxArray domain_box_array(domain_box);

        // Set the number of subdivisions for each box
        amrex::IntVect max_grid_size(params.input_max_grid_size_x, params.input_max_grid_size_y, params.input_max_grid_size_z);
        domain_box_array.maxSize(max_grid_size); // Set the maximum size of each box

        amrex::DistributionMapping distribution_mapping(domain_box_array);

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    evolve();
    amrex::Finalize();
}