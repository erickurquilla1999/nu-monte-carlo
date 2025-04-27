#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "parameters.h"
#include "particle_container.h"
#include "matter.h"
#include "moveparticles.h"

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

    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(params.input_physical_domain_size_x_cm, params.input_physical_domain_size_y_cm, params.input_physical_domain_size_z_cm)});

    amrex::Geometry geom(domain_box, &real_box);

    // Initialize the bakground matter conditions
    init_matter(matter_mfab);

    // Create the particle container
    MCParticleContainer particles(geom, distribution_mapping, domain_box_array);

    amrex::Real time_phys_s = 0.0; // speed of light in cm/s
    std::string plotfile_name;

    // Loop over the number of steps
    for (int i_step = 0; i_step < params.n_steps; ++i_step) {
        amrex::Print() << "Step: " << i_step << "\n";
        particles.Redistribute();
        particles.InsertParticles(params.test_1_n_particles, time_phys_s);
        particles.UpdateCellIndex();
        MoveParticlesMC(particles, matter_mfab, geom, params.time_step_s);

        // Write the plotfile
        if (i_step % params.write_grid == 0) {
            amrex::Print() << "Writing grid, step: " << i_step << "\n";
            plotfile_name = "plt" + std::to_string(i_step);
            amrex::WriteSingleLevelPlotfile(plotfile_name, matter_mfab, {"rho_g_ccm", "ye", "T_MeV", "IMFP_cm", "chemical_potential_MeV"}, geom, time_phys_s, i_step);
            if (i_step % params.write_particles == 0) {
                amrex::Print() << "Writing particles, step: " << i_step << "\n";
                particles.WritePlotFile(plotfile_name, "particles");
            }
        }

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