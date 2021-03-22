#ifndef __FORWARD_EULER_CPP
#define __FORWARD_EULER_CPP

#include "forward_euler.h"

namespace HyperFlow {

/* Constructor */
ForwardEulerODESolver::ForwardEulerODESolver() {};

/* Destructor */
ForwardEulerODESolver::~ForwardEulerODESolver() {};

/* Carry out one step of the ODE solution */
Vec2D ForwardEulerODESolver::operator() (const std::function<Vec2D()>& func,
                                         const Vec2D& val_to_step,
                                         const double dt)
{
    Vec2D cell_update = func();

    unsigned int cells = cell_update.size();

    for (unsigned int cell_idx=0; cell_idx<cells; cell_idx++) {
        cell_update[cell_idx] *= dt;
        cell_update[cell_idx] += val_to_step[cell_idx];
    }

    return cell_update;
}

}

#endif
