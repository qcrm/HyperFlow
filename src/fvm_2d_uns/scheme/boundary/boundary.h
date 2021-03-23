#ifndef __BOUNDARY_H
#define __BOUNDARY_H

#include <memory>

#include "../../../share/tensor/tensor.h"
#include "../../mesh/cell/cell.h"

namespace HyperFlow {

enum BoundaryCondition { Reflective = -1, Transmissive = -2, Inlet = -3 };

class BoundaryConditions
{
    
public:

    /* Constructor */
    BoundaryConditions();

    /* Constructor with inlet state */
    BoundaryConditions(std::shared_ptr<Vec1D> _inlet_state);
    
    /* Destructor */
    virtual ~BoundaryConditions();
   
    /* Calculate the boundary condition cell state */
    Vec1D calculate_boundary_condition(const BoundaryCondition& bc,
                                       const Cell& cell);

private:
    
    /* The vector of flow values */
    std::shared_ptr<Vec1D> inlet_state;

    /* Calculate the flow values for the reflective boundary condition */
    Vec1D calculate_reflective_boundary_condition(const Cell& cell);

    /* Calculate the flow values for the transmissive boundary condition */
    Vec1D calculate_transmissive_boundary_condition(const Cell& cell);

    /* Calculate the flow values for the inlet boundary condition */
    Vec1D calculate_inlet_boundary_condition(const Cell& cell);

};

}

#endif
