#ifndef __CELL_CPP
#define __CELL_CPP

#include "cell.h"

namespace HyperFlow {

/* Constructor */
Cell::Cell()
{}

/* Constructor with vertices */
Cell::Cell(const Vec2D& _vertices)
:
    Polygon(_vertices)
{}

/* Constructor with vertices and flow values */
Cell::Cell(const Vec2D& _vertices,
           const Vec1D& _flow_values)
:
    Polygon(_vertices),
    flow_values(_flow_values)
{}

/* Destructor */
Cell::~Cell()
{}

/* Get the flow values vector */
const Vec1D& Cell::get_flow_values() const
{
    return flow_values;
}

/* Set the flow values vector */
void Cell::set_flow_values(const Vec1D& _flow_values)
{
    flow_values = _flow_values;
}

}

#endif
