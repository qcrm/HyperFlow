#ifndef __BOUNDARY_CPP
#define __BOUNDARY_CPP

#include "boundary.h"

namespace HyperFlow {

/* Constructor */
BoundaryConditions::BoundaryConditions()
{}

/* Constructor with inlet state */
BoundaryConditions::BoundaryConditions(std::shared_ptr<Vec1D> _inlet_state = nullptr)
:
    inlet_state(_inlet_state)
{}

/* Destructor */
BoundaryConditions::~BoundaryConditions()
{}

/* Calculate the boundary condition cell state */
Vec1D BoundaryConditions::calculate_boundary_condition(const BoundaryCondition& bc,
                                                       const Cell& cell)
{
    Vec1D bc_flow_vals;
    if (bc == BoundaryCondition::Reflective) {
        bc_flow_vals = calculate_reflective_boundary_condition(cell);
    } else if (bc == BoundaryCondition::Transmissive) {
        bc_flow_vals = calculate_transmissive_boundary_condition(cell);
    } else if (bc == BoundaryCondition::Inlet) {
        bc_flow_vals = calculate_inlet_boundary_condition(cell);
    } else {
        std::cout << "Boundary condition '" << bc << "' not supported. Exiting." << std::endl; 
    }
    return bc_flow_vals;
}

/* Calculate the flow values for the reflective boundary condition */
Vec1D BoundaryConditions::calculate_reflective_boundary_condition(const Cell& cell)
{
    Vec1D new_flow_vals(cell.get_flow_values());
    return new_flow_vals;
}

/* Calculate the flow values for the transmissive boundary condition */
Vec1D BoundaryConditions::calculate_transmissive_boundary_condition(const Cell& cell)
{
    Vec1D new_flow_vals(cell.get_flow_values());
    return new_flow_vals;
}

/* Calculate the flow values for the inlet boundary condition */
Vec1D BoundaryConditions::calculate_inlet_boundary_condition(const Cell& cell)
{
    Vec1D new_flow_vals(*inlet_state);
    return new_flow_vals;
}

}

#endif
