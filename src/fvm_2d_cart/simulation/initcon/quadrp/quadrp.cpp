#ifndef __QUAD_RIEMANN_PROBLEM_CPP
#define __QUAD_RIEMANN_PROBLEM_CPP

#include "quadrp.h"

namespace HyperFlow {

/* Constructor */
QuadRiemannProblemInitialCondition::QuadRiemannProblemInitialCondition()
{}

/* Constructor with spatio-temporal extent, interface
 * location and all four flow states */
QuadRiemannProblemInitialCondition::QuadRiemannProblemInitialCondition(
    const double _x_left,
    const double _x_right,
    const double _y_bottom,
    const double _y_top,
    const Vec1D& _xy_interface,
    const Vec1D& _cons_nw_state,
    const Vec1D& _cons_ne_state,
    const Vec1D& _cons_sw_state,
    const Vec1D& _cons_se_state)
:
    InitialCondition(_x_left, 
                     _x_right,
                     _y_bottom,
                     _y_top,
                     _cons_nw_state.size()),
    xy_interface(_xy_interface),
    cons_nw_state(_cons_nw_state),
    cons_ne_state(_cons_ne_state),
    cons_sw_state(_cons_sw_state),
    cons_se_state(_cons_se_state)
{}
 
/* Destructor */
QuadRiemannProblemInitialCondition::~QuadRiemannProblemInitialCondition()
{}

/* Obtain the x-y coordinates of the quadrant interface */
Vec1D QuadRiemannProblemInitialCondition::get_xy_interface()
{
    return xy_interface;
}

/* Return the function for a particular dimension that
 * determines the flow values for each quadrant */
std::function<double (const double, const double)> QuadRiemannProblemInitialCondition::field_func_dim(
    const unsigned int dim
)
{
    return [this, dim](const double x, const double y) { 
        if ((x <= xy_interface[0]) && (y <= xy_interface[1])) {
            return cons_sw_state[dim];
        } else if ((x <= xy_interface[0]) && (y > xy_interface[1])) {
            return cons_nw_state[dim];
        } else if ((x > xy_interface[0]) && (y <= xy_interface[1])) {
            return cons_se_state[dim];
        } else {
            return cons_ne_state[dim];
        }
    };
}

}

#endif