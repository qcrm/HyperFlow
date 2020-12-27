#ifndef __TIMESTEP_CPP
#define __TIMESTEP_CPP

#include "timestep.h"

namespace HyperFlow {

/* Constructor */
TimeStep::TimeStep()
{}
        
/* Construct the time step with the supplied CFL
 * and pointer to the model equations */
TimeStep::TimeStep(const std::shared_ptr<Model>& _model,
                   const double _cfl)
:
    model(_model),
    cfl(_cfl)
{
    steps = 1;
}
       
/* Destructor */
TimeStep::~TimeStep()
{}

/* Calculate the timestep from the provided meshblock
 * and the CFL condition */
double TimeStep::operator() (const std::shared_ptr<Mesh>& mesh)
{
    double final_cfl = 0.0;
    double dtx = 1e16;
    double dty = 1e16;

    std::vector<std::shared_ptr<MeshBlock> > mesh_blocks = mesh->get_mesh_blocks();

    for (unsigned int b=0; b<mesh_blocks.size(); b++) {
        double dx = mesh_blocks[b]->get_dx();
        double dy = mesh_blocks[b]->get_dy();

        #pragma omp parallel for collapse(2)
        for (unsigned int i=0; i<mesh_blocks[b]->get_total_x_cells(); i++) {
            for (unsigned int j=0; j<mesh_blocks[b]->get_total_y_cells(); j++) {
                Vec1D flow_vals = mesh_blocks[b]->get_flow_values(i, j);
                Vec1D prim = model->cons_to_prim(flow_vals);
                
                double u = prim[1];
                double v = prim[2];
                double c = model->prim_speed_of_sound(prim);

                dtx = std::min(dtx, dx / (fabs(u) + c));
                dty = std::min(dty, dy / (fabs(v) + c));
            }
        }
    }

    /* This is suggested by E. Toro in order to stabilise the
     * first few iterations, particularly for challenging simulations
     * such as Riemann shock tube problems */
    if (steps <= TimeStep::initial_cfl_steps) {
        final_cfl = TimeStep::initial_cfl_multiplier * cfl;
    } else {
        final_cfl = cfl;
    }
    steps++;

    return final_cfl * std::min(dtx, dty);
}

/* Initial CFL multiplier to stabilise
 * first few iterations */
double TimeStep::initial_cfl_multiplier = 0.2;

/* Number of steps to restrict CFL for */
unsigned int TimeStep::initial_cfl_steps = 10;

}

#endif
