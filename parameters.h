#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <AMReX.H>
#include <AMReX_ParmParse.H>

struct Parameters_struct
{
    int input_n_cell_x; // Number of cells in the x-direction
    int input_n_cell_y; // Number of cells in the y-direction
    int input_n_cell_z; // Number of cells in the z-direction
    int input_max_grid_size_x; // Maximum grid size subdivisions in x of the domain for parallelization
    int input_max_grid_size_y; // Maximum grid size subdivisions in y of the domain for parallelization
    int input_max_grid_size_z; // Maximum grid size subdivisions in z of the domain for parallelization
    amrex::Real input_physical_domain_size_x_cm; // Physical size of the domain in the x-direction
    amrex::Real input_physical_domain_size_y_cm; // Physical size of the domain in the y-direction
    amrex::Real input_physical_domain_size_z_cm; // Physical size of the domain in the z-direction
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
    pp.get("input_physical_domain_size_x_cm", params.input_physical_domain_size_x_cm); // Physical size of the domain in the x-direction
    pp.get("input_physical_domain_size_y_cm", params.input_physical_domain_size_y_cm); // Physical size of the domain in the y-direction
    pp.get("input_physical_domain_size_z_cm", params.input_physical_domain_size_z_cm); // Physical size of the domain in the z-direction
}

#endif