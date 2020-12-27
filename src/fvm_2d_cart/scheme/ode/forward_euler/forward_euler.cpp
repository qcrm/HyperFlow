#ifndef __FORWARD_EULER_CPP
#define __FORWARD_EULER_CPP

#include "forward_euler.h"

namespace HyperFlow {

/* Constructor */
ForwardEulerODESolver::ForwardEulerODESolver() {};

/* Destructor */
ForwardEulerODESolver::~ForwardEulerODESolver() {};

/* Carry out one step of the ODE solution */
Vec3D ForwardEulerODESolver::operator() (const std::function<Vec3D()>& func,
                                         const Vec3D& val_to_step,
                                         const double dt)
{
    Vec3D ret_val = func();

    unsigned int x_cells = ret_val.size();
    unsigned int y_cells = ret_val[0].size();

    #pragma omp parallel for collapse(2)
    for (unsigned int i=0; i<x_cells; i++) {
        for (unsigned int j=0; j<y_cells; j++) {
            ret_val[i][j] *= dt;
            ret_val[i][j] += val_to_step[i][j];
        }
    }

    return ret_val;
}

}

#endif
