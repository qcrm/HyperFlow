#ifndef __TRANSMISSIVE_BOUNDARY_CPP
#define __TRANSMISSIVE_BOUNDARY_CPP

#include "transmissive_boundary.h"

namespace HyperFlow {

/* Constructor */
TransmissiveBoundaryCondition::TransmissiveBoundaryCondition()
{}

/* Mesh block constructor */
TransmissiveBoundaryCondition::TransmissiveBoundaryCondition(
    const std::shared_ptr<MeshBlock>& _mesh_block
)
:
    mesh_block(_mesh_block)
{}

/* Destructor */
TransmissiveBoundaryCondition::~TransmissiveBoundaryCondition()
{}

/* Apply the boundary condition to the left extent */
void TransmissiveBoundaryCondition::apply_left()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();

    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(0, j) = mesh_block->operator()(1, j);
        }
     } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(1, j) = mesh_block->operator()(2, j);
            mesh_block->operator()(0, j) = mesh_block->operator()(3, j);
        }
    }
}

// Apply the right boundary condition
void TransmissiveBoundaryCondition::apply_right()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_x_cells = mesh_block->get_total_x_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();
     
    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(total_x_cells - 1, j) = mesh_block->operator()(total_x_cells - 2, j);
        }    
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(total_x_cells - 2, j) = mesh_block->operator()(total_x_cells - 3, j);
            mesh_block->operator()(total_x_cells - 1, j) = mesh_block->operator()(total_x_cells - 4, j);
        }
    }
}

// Apply the bottom boundary codition
void TransmissiveBoundaryCondition::apply_bottom()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_x_cells = mesh_block->get_total_x_cells();

    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, 0) = mesh_block->operator()(i, 1);
        }
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, 1) = mesh_block->operator()(i, 2);
            mesh_block->operator()(i, 0) = mesh_block->operator()(i, 3);
        }
    }
}

// Apply the top boundary condition
void TransmissiveBoundaryCondition::apply_top()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_x_cells = mesh_block->get_total_x_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();
     
    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, total_y_cells - 1) = mesh_block->operator()(i, total_y_cells - 2);
        }
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, total_y_cells - 2) = mesh_block->operator()(i, total_y_cells - 3);
            mesh_block->operator()(i, total_y_cells - 1) = mesh_block->operator()(i, total_y_cells - 4);
        }
    }
}

}

#endif
