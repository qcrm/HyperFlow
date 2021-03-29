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
    double dtr = 1e16;

    CellVec1D cells = mesh->get_cells();

    unsigned int cell_size = mesh->get_cells().size();

    for (unsigned int cell_idx=0; cell_idx<cell_size; cell_idx++) {
        double dr = cells[cell_idx].radius_incircle();

        Vec1D flow_vals = cells[cell_idx].get_flow_values();
        Vec1D prim = model->cons_to_prim(flow_vals);

        double u = prim[1];
        double v = prim[2];
        double c = model->prim_speed_of_sound(prim);

        dtr = std::min(dtr, dr / (fabs(sqrt(u*u + v*v)) + c));
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

    return final_cfl * dtr;
}

/* Initial CFL multiplier to stabilise
 * first few iterations */
double TimeStep::initial_cfl_multiplier = 0.2;

/* Number of steps to restrict CFL for */
unsigned int TimeStep::initial_cfl_steps = 10;

}

#endif
