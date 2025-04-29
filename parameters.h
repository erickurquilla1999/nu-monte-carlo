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
    amrex::Real time_step_s; // Time step in seconds
    int n_steps; // Number of steps to run
    int write_grid; // Write grid every n steps
    int write_particles; // Write particles every n steps
    int n_mc_particles; // Number of Monte Carlo particles
    amrex::Real nu_Energy_center_MeV; // Energy of the neutrino in MeV
    amrex::Real nu_Energy_bottom_MeV; // Width of the neutrino energy distribution in MeV
    amrex::Real nu_Energy_top_MeV; // Width of the neutrino energy distribution in MeV
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
    pp.get("time_step_s", params.time_step_s); // Time step in seconds
    pp.get("n_steps", params.n_steps); // Number of steps to run
    pp.get("write_grid", params.write_grid); // Write grid every n steps
    pp.get("write_particles", params.write_particles); // Write particles every n steps
    pp.get("n_mc_particles", params.n_mc_particles); // Number of Monte Carlo particles
    pp.get("nu_Energy_center_MeV", params.nu_Energy_center_MeV); // Energy of the neutrino in MeV
    pp.get("nu_Energy_bottom_MeV", params.nu_Energy_bottom_MeV); // Width of the neutrino energy distribution in MeV
    pp.get("nu_Energy_top_MeV", params.nu_Energy_top_MeV); // Width of the neutrino energy distribution in MeV
}

void print_parameters(const Parameters_struct& params)
{
    amrex::Print() << "Input Parameters:" << std::endl;
    amrex::Print() << "  input_n_cell_x: " << params.input_n_cell_x << std::endl;
    amrex::Print() << "  input_n_cell_y: " << params.input_n_cell_y << std::endl;
    amrex::Print() << "  input_n_cell_z: " << params.input_n_cell_z << std::endl;
    amrex::Print() << "  input_max_grid_size_x: " << params.input_max_grid_size_x << std::endl;
    amrex::Print() << "  input_max_grid_size_y: " << params.input_max_grid_size_y << std::endl;
    amrex::Print() << "  input_max_grid_size_z: " << params.input_max_grid_size_z << std::endl;
    amrex::Print() << "  input_physical_domain_size_x_cm: " << params.input_physical_domain_size_x_cm << std::endl;
    amrex::Print() << "  input_physical_domain_size_y_cm: " << params.input_physical_domain_size_y_cm << std::endl;
    amrex::Print() << "  input_physical_domain_size_z_cm: " << params.input_physical_domain_size_z_cm << std::endl;
    amrex::Print() << "  time_step_s: " << params.time_step_s << std::endl;
    amrex::Print() << "  n_steps: " << params.n_steps << std::endl;
    amrex::Print() << "  write_grid: " << params.write_grid << std::endl;
    amrex::Print() << "  write_particles: " << params.write_particles << std::endl;
    amrex::Print() << "  n_mc_particles: " << params.n_mc_particles << std::endl;
    amrex::Print() << "  nu_Energy_center_MeV: " << params.nu_Energy_center_MeV << std::endl;
    amrex::Print() << "  nu_Energy_bottom_MeV: " << params.nu_Energy_bottom_MeV << std::endl;
    amrex::Print() << "  nu_Energy_top_MeV: " << params.nu_Energy_top_MeV << std::endl;
}

#endif