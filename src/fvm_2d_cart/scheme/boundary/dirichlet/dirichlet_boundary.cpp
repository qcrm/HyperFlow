#ifndef __DIRICHLET_BOUNDARY_CPP
#define __DIRICHLET_BOUNDARY_CPP

#include "dirichlet_boundary.h"

namespace HyperFlow {

/* Constructor */
DirichletBoundaryCondition::DirichletBoundaryCondition()
{}

/* Constructor with specified conservative flow values */
DirichletBoundaryCondition::DirichletBoundaryCondition(const std::shared_ptr<MeshBlock>& _mesh_block,
                                                       const Vec1D& _cons_flow_values)
:
    mesh_block(_mesh_block),
    cons_flow_values(_cons_flow_values)
{} 

/* Destructor */
DirichletBoundaryCondition::~DirichletBoundaryCondition()
{}

/* Apply the boundary condition to the left extent */
void DirichletBoundaryCondition::apply_left()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();

    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(0, j) = cons_flow_values;
        }
     } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(1, j) = cons_flow_values;
            mesh_block->operator()(0, j) = cons_flow_values;
        }
    }
}

// Apply the right boundary condition
void DirichletBoundaryCondition::apply_right()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_x_cells = mesh_block->get_total_x_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();
     
    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(total_x_cells - 1, j) = cons_flow_values;
        }    
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(total_x_cells - 2, j) = cons_flow_values;
            mesh_block->operator()(total_x_cells - 1, j) = cons_flow_values;
        }
    }
}

// Apply the bottom boundary codition
void DirichletBoundaryCondition::apply_bottom()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_x_cells = mesh_block->get_total_x_cells();

    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, 0) = cons_flow_values;
        }
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, 1) = cons_flow_values;
            mesh_block->operator()(i, 0) = cons_flow_values;
        }
    }
}

// Apply the top boundary condition
void DirichletBoundaryCondition::apply_top()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_x_cells = mesh_block->get_total_x_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();
     
    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, total_y_cells - 1) = cons_flow_values;
        }
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, total_y_cells - 2) = cons_flow_values;
            mesh_block->operator()(i, total_y_cells - 1) = cons_flow_values;
        }
    }
}

}

#endif
