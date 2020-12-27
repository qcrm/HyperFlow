#ifndef __GODUNOV_CPP
#define __GODUNOV_CPP

#include "godunov.h"

namespace HyperFlow {

/* Constructor */
GodunovScheme::GodunovScheme()
{}

/* Constructor with model equations, Riemann solver, time-step */
GodunovScheme::GodunovScheme(const std::shared_ptr<Model>& _model,
                             const std::shared_ptr<RiemannSolver>& _riemann,
                             const std::shared_ptr<TimeStep>& _timestep)
:
    model(_model),
    riemann(_riemann),
    timestep(_timestep)
{}

/* Destructor */
GodunovScheme::~GodunovScheme()
{}

/* Carry out the Godunov scheme on the provided mesh */
Vec3D GodunovScheme::operator() (const std::shared_ptr<MeshBlock>& mesh_block)
{
    Vec3D x_fluxes = calculate_fluxes_x(mesh_block);
    Vec3D y_fluxes = calculate_fluxes_y(mesh_block);
    Vec3D spatial_update = calculate_spatial_update(mesh_block, x_fluxes, y_fluxes);
    return spatial_update; 
}

/* Calculate the time step via the CFL condition */
double GodunovScheme::calculate_time_step(const std::shared_ptr<Mesh>& mesh)
{
    return timestep->operator()(mesh);
}

/* Number of ghost cells associated with the scheme */
const unsigned int GodunovScheme::ghost_cells = 1;

/* Calculate the fluxes in the x-direction */
Vec3D GodunovScheme::calculate_fluxes_x(const std::shared_ptr<MeshBlock>& mesh_block)
{
    unsigned int x_cells = mesh_block->get_x_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();
    unsigned int dimension = mesh_block->get_dimension();
  
    Vec3D x_fluxes(x_cells + 1, Vec2D(total_y_cells, Vec1D(dimension, 0.0)));

    for (unsigned int i=0; i<x_cells + 1; i++) {
        for (unsigned int j=0; j<total_y_cells; j++) {
            Vec1D rs_flux = riemann->operator()(
                mesh_block->get_flow_values(i, j),
                mesh_block->get_flow_values(i+1, j),
                Direction::x
            );

            x_fluxes[i][j] = rs_flux;
        }
    }
   
    return x_fluxes;
}

/* Calculate the fluxes in the y-direction */
Vec3D GodunovScheme::calculate_fluxes_y(const std::shared_ptr<MeshBlock>& mesh_block)
{
    unsigned int total_x_cells = mesh_block->get_total_x_cells();
    unsigned int y_cells = mesh_block->get_y_cells();
    unsigned int dimension = mesh_block->get_dimension();
   
    Vec3D y_fluxes(total_x_cells, Vec2D(y_cells + 1, Vec1D(dimension, 0.0)));
    
    for (unsigned int i=0; i<total_x_cells; i++) {
        for (unsigned int j=0; j<y_cells + 1; j++) {
            Vec1D rs_flux = riemann->operator()(
                mesh_block->get_flow_values(i, j),
                mesh_block->get_flow_values(i, j+1),
                Direction::y
            );

            y_fluxes[i][j] = rs_flux;
        }
    }
   
    return y_fluxes;
}

/* Calculate the spatial update via the fluxes */
Vec3D GodunovScheme::calculate_spatial_update(const std::shared_ptr<MeshBlock>& mesh_block,
                                              const Vec3D& x_fluxes,
                                              const Vec3D& y_fluxes)
{
    double one_dx = 1.0 / mesh_block->get_dx();
    double one_dy = 1.0 / mesh_block->get_dy();

    unsigned int total_x_cells = mesh_block->get_total_x_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();
    unsigned int dimension = mesh_block->get_dimension();

    Vec3D fields_update(total_x_cells, Vec2D(total_y_cells, Vec1D(dimension, 0.0)));

    for (unsigned int i=1; i<total_x_cells - 1; i++) {
        for (unsigned int j=1; j<total_y_cells - 1; j++) {
            Vec1D cell_update = (
                (one_dx * (x_fluxes[i-1][j] - x_fluxes[i][j])) + 
                (one_dy * (y_fluxes[i][j-1] - y_fluxes[i][j]))
            );

            fields_update[i][j] = cell_update;
        }
    }

    return fields_update;
}

}

#endif
