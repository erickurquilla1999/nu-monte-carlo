#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "parameters.h"
#include "particle_container.h"
#include "matter.h"
#include "particle_interactions.h"
#include "constant.h"

using namespace amrex;

void evolve()
{
    #ifdef DEBUG
        amrex::Print() << "AMReX version " << amrex::Version() << "\n";
    #endif

    // Get the simulation parameters
    Parameters_struct params;
    get_parameters(params);
    print_parameters(params);

    // Set up the simulation domain
    amrex::IntVect domain_low_cells(0,0,0);
    amrex::IntVect domain_high_cells(params.input_n_cell_x-1, params.input_n_cell_y-1, params.input_n_cell_z-1);

    // Create the domain box
    amrex::Box domain_box(domain_low_cells, domain_high_cells);
    amrex::BoxArray domain_box_array(domain_box);

    // Set the number of subdivisions for each box
    amrex::IntVect max_grid_size(params.input_max_grid_size_x, params.input_max_grid_size_y, params.input_max_grid_size_z);
    domain_box_array.maxSize(max_grid_size); // Set the maximum size of each box

    // Create the distribution mapping
    amrex::DistributionMapping distribution_mapping(domain_box_array);

    const int num_boxes = domain_box_array.size();

    #ifdef DEBUG
        amrex::Print() << "Number of boxes: " << num_boxes << "\n";
    #endif

    amrex::MultiFab matter_mfab(domain_box_array, distribution_mapping, MatterData::ncomps, 0);
    amrex::MultiFab nu_mfab(domain_box_array, distribution_mapping, NeutrinoData::ncomps, 0);

    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(params.input_physical_domain_size_x_cm, params.input_physical_domain_size_y_cm, params.input_physical_domain_size_z_cm)});

    amrex::Geometry geom(domain_box, &real_box);

    const amrex::Real* dx = geom.CellSize();
    const amrex::Real* plo = geom.ProbLo();
    amrex::Real dx_local[3] = {dx[0], dx[1], dx[2]};
    amrex::Real plo_local[3] = {plo[0], plo[1], plo[2]};
    amrex::Real cellvolume = dx[0] * dx[1] * dx[2];

    // Initialize the bakground matter conditions
    init_matter(matter_mfab, geom, params.simulation_type);

    // Create the particle container
    MCParticleContainer particles(geom, distribution_mapping, domain_box_array);

    amrex::Real time_phys_s = 0.0; // speed of light in cm/s
    std::string plotfile_name;

    const amrex::Real dtdE3_3dOmegadx3 = ( params.time_step_s ) * //dt
                                ( (params.nu_Energy_top_MeV*params.nu_Energy_top_MeV*params.nu_Energy_top_MeV - params.nu_Energy_bottom_MeV*params.nu_Energy_bottom_MeV*params.nu_Energy_bottom_MeV)/3.0 ) * //dE3
                                ( 4.0 * MathConst::pi ) * // dOmega
                                ( cellvolume ); // dx3

    amrex::Real n_nu_per_mc_particles = 0.0;
    compute_nu_per_MC_particles(matter_mfab, params.n_mc_particles, n_nu_per_mc_particles, params.nu_Energy_center_MeV, dtdE3_3dOmegadx3);
    amrex::Print() << "Number of neutrinos per MC particle: " << n_nu_per_mc_particles << "\n";
    amrex::Print() << "dtdE3_3dOmegadx3: " << dtdE3_3dOmegadx3 << "\n";

    // Write matter plotfile
    amrex::WriteSingleLevelPlotfile("pltmatter", matter_mfab, {"rho_g_ccm", "ye", "T_MeV", "IMFP_cm", "chemical_potential_MeV"}, geom, 0.0, 0);

    // Loop over the number of steps
    for (int i_step = 0; i_step < params.n_steps+1; ++i_step) {
        amrex::Print() << "Step: " << i_step << "\n";
        particles.EmissionParticles(matter_mfab, n_nu_per_mc_particles, params.nu_Energy_center_MeV, dtdE3_3dOmegadx3, time_phys_s, params.simulation_type);

        // Write the plotfile
        if (i_step % params.write_grid == 0) {
            amrex::Print() << "Writing grid, step: " << i_step << "\n";
            compute_nu_n_and_f(particles, geom, nu_mfab);
            plotfile_name = "plt" + std::to_string(i_step);
            amrex::WriteSingleLevelPlotfile(plotfile_name, nu_mfab, {"n_invcm3", "fx_invcm3", "fy_invcm3", "fz_invcm3"}, geom, time_phys_s, i_step);
            if (i_step % params.write_particles == 0) {
                amrex::Print() << "Writing particles, step: " << i_step << "\n";
                particles.WritePlotFile(plotfile_name, "particles");
            }
        }
        
        MoveParticlesMC(particles, matter_mfab, geom, params.time_step_s);
        particles.AbsorptionParticles();

        // Update the time
        time_phys_s += params.time_step_s;

    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    evolve();
    amrex::Finalize();
}