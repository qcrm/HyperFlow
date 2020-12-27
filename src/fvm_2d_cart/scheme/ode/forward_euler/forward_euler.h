#ifndef __FORWARD_EULER_H
#define __FORWARD_EULER_H

#include <omp.h>

#include "../ode/ode.h"

namespace HyperFlow {

class ForwardEulerODESolver
:
    public ODESolver
{
    
public:
    
    /* Constructor */
    ForwardEulerODESolver();

    /* Destructor */
    virtual ~ForwardEulerODESolver();
    
    /* Carry out one step of the ODE solution */
    virtual Vec3D operator() (const std::function<Vec3D()>& func,
                              const Vec3D& val_to_step,
                              const double dt);

};

}

#endif
