#ifndef __FORWARD_EULER_H
#define __FORWARD_EULER_H

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
    virtual Vec2D operator() (const std::function<Vec2D()>& func,
                              const Vec2D& val_to_step,
                              const double dt);

};

}

#endif
