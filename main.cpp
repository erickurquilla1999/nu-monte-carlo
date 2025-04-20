#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "parameters.h"
#include "particle_container.h"
#include "matter.h"

using namespace amrex;

void evolve()
{
    amrex::Print() << "AMReX version " << amrex::Version() << "\n";

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
    amrex::Print() << "Number of boxes: " << num_boxes << "\n";

    amrex::MultiFab matter_mfab(domain_box_array, distribution_mapping, MatterData::ncomps, 0);

    amrex::RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                            {AMREX_D_DECL(params.input_physical_domain_size_x_cm, params.input_physical_domain_size_y_cm, params.input_physical_domain_size_z_cm)});

    amrex::Geometry geom(domain_box, &real_box);
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> cell_size_cm = geom.CellSizeArray();

    // Initialize the bakground matter conditions
    init_matter(matter_mfab);

    // Create the particle container
    MCParticleContainer particles(geom, distribution_mapping, domain_box_array);
    particles.InsertParticles(params.test_1_n_particles);

    // print particles value
    particles.LoopParticlesPrint();

    // Loop over the number of steps
    for (int i_step = 0; i_step < params.n_steps; ++i_step) {
        amrex::Print() << "Step: " << i_step << "\n";
        particles.MoveParticles(params.time_step_s);
        particles.InsertParticles(params.test_1_n_particles);
        particles.Redistribute();
        particles.LoopParticlesPrint();
    }

    // WriteSingleLevelPlotfile("plt001", matter_mfab, {"matter_mfab"}, geom, 0.0, 0);
    // particles.Checkpoint("plt001", "particle0");

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    evolve();
    amrex::Finalize();
}