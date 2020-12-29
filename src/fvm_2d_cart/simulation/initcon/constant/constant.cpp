#ifndef __CONSTANT_INITCON_CPP
#define __CONSTANT_INITCON_CPP

#include "constant.h"

namespace HyperFlow {

/* Constructor */
ConstantInitialCondition::ConstantInitialCondition()
{}

/* Constructor with constant initial conservative flow state */
ConstantInitialCondition::ConstantInitialCondition(
    const double _x_left,
    const double _x_right,
    const double _y_bottom,
    const double _y_top,
    const Vec1D& _cons_init_state)
:
    InitialCondition(_x_left, 
                     _x_right,
                     _y_bottom,
                     _y_top,
                     _cons_init_state.size()),
    cons_init_state(_cons_init_state)
{}
 
/* Destructor */
ConstantInitialCondition::~ConstantInitialCondition()
{}

/* Return the function for a particular dimension that
         * determines the constant flow values */
std::function<double (const double, const double)> ConstantInitialCondition::field_func_dim(
    const unsigned int dim
) {
    return [this, dim](const double x, const double y) { 
        return cons_init_state[dim];
    };
}

}

#endif