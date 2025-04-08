#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

struct Parameters_struct
{
    int input_n_cell_x; // Number of cells in the x-direction
    int input_n_cell_y; // Number of cells in the y-direction
    int input_n_cell_z; // Number of cells in the z-direction
    int input_max_grid_size_x; // Maximum grid size subdivisions in x of the domain for parallelization
    int input_max_grid_size_y; // Maximum grid size subdivisions in y of the domain for parallelization
    int input_max_grid_size_z; // Maximum grid size subdivisions in z of the domain for parallelization
};

void get_parameters(Parameters_struct& params)
{
    amrex::ParmParse pp;
    pp.get("input_n_cell_x", params.input_n_cell_x); // Number of cells in the x-direction
    pp.get("input_n_cell_y", params.input_n_cell_y); // Number of cells in the y-direction
    pp.get("input_n_cell_z", params.input_n_cell_z); // Number of cells in the z-direction
    pp.get("input_max_grid_size_x", params.input_max_grid_size_x); // Maximum grid size subdivisions in x of the domain for parallelization
    pp.get("input_max_grid_size_y", params.input_max_grid_size_y); // Maximum grid size subdivisions in y of the domain for parallelization 
    pp.get("input_max_grid_size_z", params.input_max_grid_size_z); // Maximum grid size subdivisions in z of the domain for parallelization
}

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