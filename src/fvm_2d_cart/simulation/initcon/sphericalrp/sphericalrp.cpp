#ifndef __SPHERICAL_RIEMANN_PROBLEM_CPP
#define __SPHERICAL_RIEMANN_PROBLEM_CPP

#include "sphericalrp.h"

namespace HyperFlow {

/* Constructor */
SphericalRiemannProblemInitialCondition::SphericalRiemannProblemInitialCondition()
{}

/* Constructor with spatio-temporal extent, sphere
 * geometry and internal and external flow states */
SphericalRiemannProblemInitialCondition::SphericalRiemannProblemInitialCondition(
    const double _x_left,
    const double _x_right,
    const double _y_bottom,
    const double _y_top,
    const double _radius,
    const double _x_origin,
    const double _y_origin,
    const Vec1D& _cons_sphere_state,
    const Vec1D& _cons_ext_state)
:
    InitialCondition(_x_left, 
                     _x_right,
                     _y_bottom,
                     _y_top,
                     _cons_sphere_state.size()),
    radius(_radius),
    x_origin(_x_origin),
    y_origin(_y_origin),
    cons_sphere_state(_cons_sphere_state),
    cons_ext_state(_cons_ext_state)
{}
 
/* Destructor */
SphericalRiemannProblemInitialCondition::~SphericalRiemannProblemInitialCondition()
{}

/* Obtain the radius of the sphere flow state */
double SphericalRiemannProblemInitialCondition::get_radius()
{
    return radius;
}

/* Obtain the x coordinate of the sphere flow state centre */
double SphericalRiemannProblemInitialCondition::get_x_origin()
{
    return x_origin;
}

/* Obtain the y coordinate of the sphere flow state centre */
double SphericalRiemannProblemInitialCondition::get_y_origin()
{
    return y_origin;
}

/* Return the function for a particular dimension that
 * determines the flow values internal and external
 * to the spherical region */
std::function<double (const double, const double)> SphericalRiemannProblemInitialCondition::field_func_dim(
    const unsigned int dim
) {
    return [this, dim](const double x, const double y) { 
        if (((x - x_origin)*(x - x_origin) + (y - y_origin)*(y - y_origin)) <= (radius*radius)) {
            return cons_sphere_state[dim];
        } else {
            return cons_ext_state[dim];
        }
    };
}

}

#endif