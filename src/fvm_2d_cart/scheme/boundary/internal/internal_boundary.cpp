#ifndef __INTERNAL_BOUNDARY_CPP
#define __INTERNAL_BOUNDARY_CPP

#include "internal_boundary.h"

namespace HyperFlow {

/* Constructor */
InternalBoundaryCondition::InternalBoundaryCondition()
{}

/* Mesh block constructor */
InternalBoundaryCondition::InternalBoundaryCondition(
    const std::shared_ptr<MeshBlock>& _internal_block,
    const std::shared_ptr<MeshBlock>& _external_block
)
:
    internal_block(_internal_block),
    external_block(_external_block)
{}

/* Destructor */
InternalBoundaryCondition::~InternalBoundaryCondition()
{}

/* Apply the boundary condition to the left extent */
void InternalBoundaryCondition::apply_left()
{
    unsigned int ghost_cells = internal_block->get_ghost_cells();
    unsigned int internal_y_cells = internal_block->get_total_y_cells();
    unsigned int external_x_cells = external_block->get_total_x_cells();

    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<internal_y_cells; j++) {
            internal_block->operator()(0, j) = external_block->operator()(external_x_cells - 2, j);
        }
     } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<internal_y_cells; j++) {
            internal_block->operator()(1, j) = external_block->operator()(external_x_cells - 3, j);
            internal_block->operator()(0, j) = external_block->operator()(external_x_cells - 4, j);
        }
    }
}

// Apply the right boundary condition
void InternalBoundaryCondition::apply_right()
{
    unsigned int ghost_cells = internal_block->get_ghost_cells();
    unsigned int internal_x_cells = internal_block->get_total_x_cells();
    unsigned int internal_y_cells = internal_block->get_total_y_cells();
     
    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<internal_y_cells; j++) {
            internal_block->operator()(internal_x_cells - 1, j) = external_block->operator()(1, j);
        }    
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<internal_y_cells; j++) {
            internal_block->operator()(internal_x_cells - 1, j) = external_block->operator()(3, j);
            internal_block->operator()(internal_x_cells - 2, j) = external_block->operator()(2, j);
        }
    }
}

// Apply the bottom boundary codition
void InternalBoundaryCondition::apply_bottom()
{
    unsigned int ghost_cells = internal_block->get_ghost_cells();
    unsigned int internal_x_cells = internal_block->get_total_x_cells();
    unsigned int external_y_cells = external_block->get_total_y_cells();

    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<internal_x_cells; i++) {
            internal_block->operator()(i, 0) = external_block->operator()(i, external_y_cells - 2);
        }
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<internal_x_cells; i++) {
            internal_block->operator()(i, 1) = external_block->operator()(i, external_y_cells - 3);
            internal_block->operator()(i, 0) = external_block->operator()(i, external_y_cells - 4);
        }
    }
}

// Apply the top boundary condition
void InternalBoundaryCondition::apply_top()
{
    unsigned int ghost_cells = internal_block->get_ghost_cells();
    unsigned int internal_x_cells = internal_block->get_total_x_cells();
    unsigned int internal_y_cells = internal_block->get_total_y_cells();
     
    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<internal_x_cells; i++) {
            internal_block->operator()(i, internal_y_cells - 1) = external_block->operator()(i, 1);
        }
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<internal_x_cells; i++) {
            internal_block->operator()(i, internal_y_cells - 1) = external_block->operator()(i, 3);
            internal_block->operator()(i, internal_y_cells - 2) = external_block->operator()(i, 2);
        }
    }
}

}

#endif
