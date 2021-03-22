#ifndef __INITCON_CPP
#define __INITCON_CPP

#include "initcon.h"

namespace HyperFlow {

/* Constructor */
InitialCondition::InitialCondition()
{}

/* Constructor from spatial and temporal extents */
InitialCondition::InitialCondition(const double _x_left,
                                   const double _x_right,
                                   const double _y_bottom,
                                   const double _y_top,
                                   const double _t_start,
                                   const double _t_end,
                                   const unsigned int _dimension)
:
    x_left(_x_left),
    x_right(_x_right),
    y_bottom(_y_bottom),
    y_top(_y_top),
    t_start(_t_start),
    t_end(_t_end),
    dimension(_dimension)
{}

/* Destructor */
InitialCondition::~InitialCondition()
{}

/* Obtain the initial flow values at the
 * provided coordinates */
Vec1D InitialCondition::operator() (const double x, const double y)
{
    Vec1D func_vals;

    for (unsigned int dim=0; dim<dimension; dim++) {
        func_vals.push_back(field_func_dim(dim)(x, y));
    }
    
    return func_vals;
}

/* Obtain the left extent of the simulation */
double InitialCondition::get_x_left()
{
    return x_left;
}

/* Obtain the right extent of the simulation */
double InitialCondition::get_x_right()
{
    return x_right;
}

/* Obtain the bottom extent of the simulation */
double InitialCondition::get_y_bottom()
{
    return y_bottom;
}

/* Obtain the top extent of the simulation */
double InitialCondition::get_y_top()
{
    return y_top;
}

/* Obtain the simualation start time */
double InitialCondition::get_t_start()
{
    return t_start;
}

/* Obtain the simulation end time */ 
double InitialCondition::get_t_end()
{
    return t_end;
}

/* Obtain the dimensionality of the model system */ 
unsigned int InitialCondition::get_dimension()
{
    return dimension;
}

}

#endif
