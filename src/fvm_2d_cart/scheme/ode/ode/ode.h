#ifndef __ODE_H
#define __ODE_H

#include <functional>

#include "../../../tensor/tensor.h"

namespace HyperFlow {

class ODESolver
{

public:

    /* Constructor */
    ODESolver();

    /* Destructor */
    virtual ~ODESolver();

    /* Carry out one step of the ODE solution */
    virtual Vec3D operator() (const std::function<Vec3D()>& func,
                              const Vec3D& val_to_step,
                              const double dt) = 0;

};

}

#endif
