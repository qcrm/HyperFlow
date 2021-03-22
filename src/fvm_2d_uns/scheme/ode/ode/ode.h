#ifndef __ODE_H
#define __ODE_H

#include <functional>

#include "../../../../share/tensor/tensor.h"

namespace HyperFlow {

class ODESolver
{

public:

    /* Constructor */
    ODESolver();

    /* Destructor */
    virtual ~ODESolver();

    /* Carry out one step of the ODE solution */
    virtual Vec2D operator() (const std::function<Vec2D()>& func,
                              const Vec2D& val_to_step,
                              const double dt) = 0;

};

}

#endif
