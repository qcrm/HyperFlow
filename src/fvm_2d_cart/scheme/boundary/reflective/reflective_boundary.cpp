#ifndef __REFLECTIVE_BOUNDARY_CPP
#define __REFLECTIVE_BOUNDARY_CPP

#include "reflective_boundary.h"

namespace HyperFlow {

/* Constructor */
ReflectiveBoundaryCondition::ReflectiveBoundaryCondition()
{}

/* Mesh block constructor */
ReflectiveBoundaryCondition::ReflectiveBoundaryCondition(
    const std::shared_ptr<MeshBlock>& _mesh_block
)
:
    mesh_block(_mesh_block)
{}

/* Destructor */
ReflectiveBoundaryCondition::~ReflectiveBoundaryCondition()
{}

/* Apply the boundary condition to the left extent */
void ReflectiveBoundaryCondition::apply_left()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();

    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(0, j) = mesh_block->operator()(1, j);
            mesh_block->operator()(0, j)[1] *= -1.0;
        }
     } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(1, j) = mesh_block->operator()(2, j);
            mesh_block->operator()(0, j) = mesh_block->operator()(3, j);

            mesh_block->operator()(1, j)[1] *= -1.0;
            mesh_block->operator()(0, j)[1] *= -1.0;
        }
    }
}

// Apply the right boundary condition
void ReflectiveBoundaryCondition::apply_right()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_x_cells = mesh_block->get_total_x_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();
     
    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(total_x_cells - 1, j) = mesh_block->operator()(total_x_cells - 2, j);
            mesh_block->operator()(total_x_cells - 1, j)[1] *= -1.0;
        }    
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int j=0; j<total_y_cells; j++) {
            mesh_block->operator()(total_x_cells - 2, j) = mesh_block->operator()(total_x_cells - 3, j);
            mesh_block->operator()(total_x_cells - 1, j) = mesh_block->operator()(total_x_cells - 4, j);

            mesh_block->operator()(total_x_cells - 2, j)[1] *= -1.0;
            mesh_block->operator()(total_x_cells - 1, j)[1] *= -1.0;
        }
    }
}

// Apply the bottom boundary codition
void ReflectiveBoundaryCondition::apply_bottom()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_x_cells = mesh_block->get_total_x_cells();

    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, 0) = mesh_block->operator()(i, 1);
            mesh_block->operator()(i, 0)[2] *= -1.0;
        }
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, 1) = mesh_block->operator()(i, 2);
            mesh_block->operator()(i, 0) = mesh_block->operator()(i, 3);

            mesh_block->operator()(i, 1)[2] *= -1.0;
            mesh_block->operator()(i, 0)[2] *= -1.0;
        }
    }
}

// Apply the top boundary condition
void ReflectiveBoundaryCondition::apply_top()
{
    unsigned int ghost_cells = mesh_block->get_ghost_cells();
    unsigned int total_x_cells = mesh_block->get_total_x_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();
     
    if (ghost_cells == 1) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, total_y_cells - 1) = mesh_block->operator()(i, total_y_cells - 2);
            mesh_block->operator()(i, total_y_cells - 1)[2] *= -1.0;
        }
    } else if (ghost_cells == 2) {
        #pragma omp parallel for collapse(1)
        for (unsigned int i=0; i<total_x_cells; i++) {
            mesh_block->operator()(i, total_y_cells - 2) = mesh_block->operator()(i, total_y_cells - 3);
            mesh_block->operator()(i, total_y_cells - 1) = mesh_block->operator()(i, total_y_cells - 4);

            mesh_block->operator()(i, total_y_cells - 2)[2] *= -1.0;
            mesh_block->operator()(i, total_y_cells - 1)[2] *= -1.0;
        }
    }
}

}

#endif
