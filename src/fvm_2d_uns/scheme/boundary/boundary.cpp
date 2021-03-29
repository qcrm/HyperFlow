#ifndef __BOUNDARY_CPP
#define __BOUNDARY_CPP

#include "boundary.h"

namespace HyperFlow {

/* Constructor */
BoundaryConditions::BoundaryConditions()
{}

/* Constructor with equation model */
BoundaryConditions::BoundaryConditions(std::shared_ptr<Model> _model)
{}

/* Constructor with equation model and inlet state */
BoundaryConditions::BoundaryConditions(std::shared_ptr<Model> _model,
                                       std::shared_ptr<Vec1D> _inlet_state = nullptr)
:
    model(_model),
    inlet_state(_inlet_state)
{}

/* Destructor */
BoundaryConditions::~BoundaryConditions()
{}

/* Calculate the boundary condition cell state */
Vec1D BoundaryConditions::calculate_boundary_condition(const BoundaryCondition& bc,
                                                       const Cell& cell,
                                                       const double edge_normal_angle)
{
    Vec1D bc_flow_vals;
    if (bc == BoundaryCondition::Reflective) {
        bc_flow_vals = calculate_reflective_boundary_condition(cell, edge_normal_angle);
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
Vec1D BoundaryConditions::calculate_reflective_boundary_condition(const Cell& cell,
                                                                  const double edge_normal_angle)
{
    Vec1D new_flow_vals(cell.get_flow_values());
    Vec1D rot_nbrc_values = model->rotate_flow_values(edge_normal_angle, new_flow_vals);
    rot_nbrc_values[1] *= -1.0;
    Vec1D rot_back_nbrc_values = model->rotate_flow_values(-1.0*edge_normal_angle, rot_nbrc_values);
    return rot_back_nbrc_values;
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
