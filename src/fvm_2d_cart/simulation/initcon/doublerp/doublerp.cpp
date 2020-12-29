#ifndef __DOUBLE_RIEMANN_PROBLEM_CPP
#define __DOUBLE_RIEMANN_PROBLEM_CPP

#include "doublerp.h"

namespace HyperFlow {

/* Constructor */
DoubleRiemannProblemInitialCondition::DoubleRiemannProblemInitialCondition()
{}

/* Parameterised constructor */
DoubleRiemannProblemInitialCondition::DoubleRiemannProblemInitialCondition(
    const double _x_left,
    const double _x_right,
    const double _y_bottom,
    const double _y_top,
    const double _x_lm_interface,
    const double _x_mr_interface,
    const Vec1D& _cons_left_state,
    const Vec1D& _cons_middle_state,
    const Vec1D& _cons_right_state)
:
    InitialCondition(_x_left, 
                     _x_right,
                     _y_bottom,
                     _y_top,
                     _cons_left_state.size()),
    x_lm_interface(_x_lm_interface),
    x_mr_interface(_x_mr_interface),
    cons_left_state(_cons_left_state),
    cons_middle_state(_cons_middle_state),
    cons_right_state(_cons_right_state)
{}
 
/* Destructor */
DoubleRiemannProblemInitialCondition::~DoubleRiemannProblemInitialCondition()
{}

/* Return the function for a particular dimension that
 * determines the flow values for each state */
std::function<double (const double, const double)> DoubleRiemannProblemInitialCondition::field_func_dim(
    const unsigned int dim
)
{
    return [this, dim](const double x, const double y) { 
        if (x <= x_lm_interface) {
            return cons_left_state[dim];
        } else if ((x > x_lm_interface) && (x <= x_mr_interface)) {
            return cons_middle_state[dim];
        } else {
            return cons_right_state[dim];
        }
    };
}

}

#endif