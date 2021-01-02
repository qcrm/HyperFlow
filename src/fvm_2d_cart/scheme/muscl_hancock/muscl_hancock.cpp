#ifndef __MUSCL_HANCOCK_CPP
#define __MUSCL_HANCOCK_CPP

#include "muscl_hancock.h"

namespace HyperFlow {

double sign(const double val)
{
    return ((val > 0) - (val < 0));    
};

/* Constructor */
MUSCLHancockScheme::MUSCLHancockScheme()
{}; 

/* Constructor with model equations, Riemann solver, time-step, limiter */
MUSCLHancockScheme::MUSCLHancockScheme(
    const std::shared_ptr<Model> _model,
    const std::shared_ptr<RiemannSolver> _riemann,
    const std::shared_ptr<TimeStep> _time_step,
    const std::shared_ptr<Limiter> _limiter,
    const std::shared_ptr<Mesh> _mesh
)
:
    model(_model),
    riemann(_riemann),
    time_step(_time_step),
    limiter(_limiter),
    mesh(_mesh)
{
    tol = 1e-6;
    omega = 0.0;
}

/* Destructor */
MUSCLHancockScheme::~MUSCLHancockScheme()
{};

/* Carry out the MUSCL-Hancock scheme on the provided mesh block */
Vec3D MUSCLHancockScheme::operator() (const std::shared_ptr<MeshBlock>& mesh_block)
{
    Vec3D slopes_x = calculate_slopes_x(mesh_block);
    Vec3D slopes_y = calculate_slopes_y(mesh_block);

    Vec4D reconstruction = calculate_boundary_extrapolations(mesh_block, slopes_x, slopes_y);
    evolve_boundary_reconstruction(mesh_block, reconstruction);

    Vec3D x_fluxes = calculate_fluxes_x(mesh_block, reconstruction);
    Vec3D y_fluxes = calculate_fluxes_y(mesh_block, reconstruction);

    Vec3D spatial_update = calculate_spatial_update(mesh_block, x_fluxes, y_fluxes);

    return spatial_update; 
}

/* Calculate the time step via the CFL condition */
double MUSCLHancockScheme::calculate_time_step(const std::shared_ptr<Mesh>& full_mesh)
{
    return time_step->operator()(full_mesh);
}

/* Number of ghost cells associated with the scheme */
const unsigned int MUSCLHancockScheme::ghost_cells = 2;

/* Calculate the x-slope limiter values for each dimension */
Vec1D MUSCLHancockScheme::calculate_slope_limiters_x(const std::shared_ptr<MeshBlock>& mesh_block,
                                                     const unsigned int x_index,
                                                     const unsigned int y_index)
{
    unsigned int dimension = mesh_block->get_dimension();
    
    Vec1D limiters(dimension, 0.0);
    Vec1D dleft = mesh_block->delta_x_back(x_index, y_index);
    Vec1D dright = mesh_block->delta_x_forward(x_index, y_index);

    for (unsigned int dim=0; dim<dimension; dim++) {
        if (fabs(dleft[dim]) <= tol) {
            if (dleft[dim] == 0.0) {
                dleft[dim]= tol;
            } else {
                dleft[dim] = tol * sign(dleft[dim]);
            }
        }

        if (fabs(dright[dim]) <= tol) {
            if (dright[dim] == 0.0) {
                dright[dim] = tol;
            } else {
                dright[dim] = tol * sign(dright[dim]);
            }
        }

        double r = dleft[dim] / dright[dim];
        limiters[dim] = limiter->operator()(r);
    }

    return limiters;
}

/* Calculate the y-slope limiter values for each dimension */
Vec1D MUSCLHancockScheme::calculate_slope_limiters_y(const std::shared_ptr<MeshBlock>& mesh_block,
                                                     const unsigned int x_index,
                                                     const unsigned int y_index)
{
    unsigned int dimension = mesh_block->get_dimension();
    
    Vec1D limiters(dimension, 0.0);
    Vec1D dbottom = mesh_block->delta_y_back(x_index, y_index);
    Vec1D dtop = mesh_block->delta_y_forward(x_index, y_index);

    for (unsigned int dim = 0; dim < dimension; dim++) {
        if (fabs(dbottom[dim]) <= tol) {
            if (dbottom[dim] == 0.0) {
                dbottom[dim]= tol;
            } else {
                dbottom[dim] = tol * sign(dbottom[dim]);
            }
        }

        if (fabs(dtop[dim]) <= tol) {
            if (dtop[dim] == 0.0) {
                dtop[dim] = tol;
            } else {
                dtop[dim] = tol * sign(dtop[dim]);
            }
        }

        double r = dbottom[dim] / dtop[dim];
        limiters[dim] = limiter->operator()(r);
    }

    return limiters;
}

/* Calculate the x-direction slopes for the linear reconstruction */
Vec3D MUSCLHancockScheme::calculate_slopes_x(const std::shared_ptr<MeshBlock>& mesh_block)
{
    unsigned int total_cells_x = mesh_block->get_total_x_cells();
    unsigned int total_cells_y = mesh_block->get_total_y_cells();
    unsigned int dimension = mesh_block->get_dimension();
    
    Vec3D slopes_x(total_cells_x, Vec2D(total_cells_y, Vec1D(dimension, 0.0)));
   
    #pragma omp parallel for collapse(2)
    for (unsigned int i=1; i<total_cells_x - 1; i++) {
        for (unsigned int j=1; j<total_cells_y - 1; j++) {
            Vec1D limiters;
            
            if (limiter) {
                limiters = calculate_slope_limiters_x(mesh_block, i, j);
            } else {
                for (unsigned int dim=0; dim<dimension; dim++) {
                    limiters.push_back(0.0);
                }
            }

            Vec1D d_slopes_left(dimension, 0.0);
            Vec1D d_slopes_right(dimension, 0.0);

            for (unsigned int dim=0; dim<dimension; dim++) {
                d_slopes_left[dim] = 0.5 * (1.0 + omega) * mesh_block->delta_x_back(i, j)[dim];
                d_slopes_right[dim] = 0.5 * (1.0 - omega) * mesh_block->delta_x_forward(i, j)[dim];
                slopes_x[i][j][dim] = limiters[dim] * (d_slopes_left[dim] + d_slopes_right[dim]);
            }
        }
    }
    
    return slopes_x;
}

/* Calculate the y-direction slopes for the linear reconstruction */
Vec3D MUSCLHancockScheme::calculate_slopes_y(const std::shared_ptr<MeshBlock>& mesh_block)
{
    unsigned int total_cells_x = mesh_block->get_total_x_cells();
    unsigned int total_cells_y = mesh_block->get_total_y_cells();   
    unsigned int dimension = mesh_block->get_dimension();
    
    Vec3D slopes_y(total_cells_x, Vec2D(total_cells_y, Vec1D(dimension, 0.0)));
   
    #pragma omp parallel for collapse(2)
    for (unsigned int i=1; i<total_cells_x - 1; i++) {
        for (unsigned int j=1; j<total_cells_y - 1; j++) {
            Vec1D limiters;

            if (limiter) {
                limiters = calculate_slope_limiters_y(mesh_block, i, j);
            } else {
                for (unsigned int dim=0; dim<dimension; dim++) {
                    limiters.push_back(0.0);
                }
            }

            Vec1D d_slopes_bottom(dimension, 0.0);
            Vec1D d_slopes_top(dimension, 0.0);

            for (unsigned int dim=0; dim<dimension; dim++) {
                d_slopes_bottom[dim] = 0.5 * (1.0 + omega) * mesh_block->delta_y_back(i, j)[dim];
                d_slopes_top[dim] = 0.5 * (1.0 - omega) * mesh_block->delta_y_forward(i, j)[dim];
                slopes_y[i][j][dim] = limiters[dim] * (d_slopes_bottom[dim] + d_slopes_top[dim]);
            }
        }
    }
    
    return slopes_y;
}

/* Calculate the linear boundary extrapolations */
Vec4D MUSCLHancockScheme::calculate_boundary_extrapolations(const std::shared_ptr<MeshBlock>& mesh_block,
                                                            const Vec3D& slopes_x,
                                                            const Vec3D& slopes_y)
{
    unsigned int total_cells_x = mesh_block->get_total_x_cells();
    unsigned int total_cells_y = mesh_block->get_total_y_cells();
    unsigned int dimension = mesh_block->get_dimension();
    
    Vec4D boundary_extrapolations(total_cells_x, Vec3D(total_cells_y, Vec2D(dimension, Vec1D(4, 0.0))));

    #pragma omp parallel for collapse(2)
    for (unsigned int i=1; i<total_cells_x - 1; i++) {
        for (unsigned int j=1; j<total_cells_y - 1; j++) {
            Vec1D cell_avgs = mesh_block->operator()(i, j);
            Vec1D cell_slopes_x = slopes_x[i][j];
            Vec1D cell_slopes_y = slopes_y[i][j];

            for (unsigned int dim=0; dim<dimension; dim++) {
                Vec1D cell_extraps(4, 0.0);

                cell_extraps[0] = cell_avgs[dim] - 0.5 * cell_slopes_x[dim];
                cell_extraps[1] = cell_avgs[dim] + 0.5 * cell_slopes_x[dim];
                cell_extraps[2] = cell_avgs[dim] - 0.5 * cell_slopes_y[dim];
                cell_extraps[3] = cell_avgs[dim] + 0.5 * cell_slopes_y[dim];

                boundary_extrapolations[i][j][dim] = cell_extraps;
            }
        }
    }

    return boundary_extrapolations;
}

/* Predictor-corrector half time-step estimate of update */
void MUSCLHancockScheme::evolve_boundary_reconstruction(const std::shared_ptr<MeshBlock>& mesh_block,
                                                        Vec4D& reconstruction)
{
    double dt = calculate_time_step(mesh);

    double half_dt_dx = 0.5 * dt / mesh_block->get_dx();
    double half_dt_dy = 0.5 * dt / mesh_block->get_dy();

    unsigned int dimension = mesh_block->get_dimension();

    #pragma omp parallel for collapse(2)
    for (unsigned int i=1; i<mesh_block->get_total_x_cells() - 1; i++) {
        for (unsigned int j=1; j<mesh_block->get_total_y_cells() - 1; j++) {
            Vec1D rec_lft(dimension, 0.0);
            Vec1D rec_rgt(dimension, 0.0);
            Vec1D rec_bot(dimension, 0.0);
            Vec1D rec_top(dimension, 0.0);

            for (unsigned int dim=0; dim<dimension; dim++) {
                rec_lft[dim] = reconstruction[i][j][dim][0];
                rec_rgt[dim] = reconstruction[i][j][dim][1];
                rec_bot[dim] = reconstruction[i][j][dim][2];
                rec_top[dim] = reconstruction[i][j][dim][3];
            }

            Vec1D flux_lft = model->cons_flux_x(rec_lft);
            Vec1D flux_rgt = model->cons_flux_x(rec_rgt);
            Vec1D flux_bot = model->cons_flux_y(rec_bot);
            Vec1D flux_top = model->cons_flux_y(rec_top);

            for (unsigned int dim = 0; dim < dimension; dim++) {
                reconstruction[i][j][dim][0] += half_dt_dx * (flux_lft[dim] - flux_rgt[dim]);
                reconstruction[i][j][dim][1] += half_dt_dx * (flux_lft[dim] - flux_rgt[dim]);
                reconstruction[i][j][dim][2] += half_dt_dy * (flux_bot[dim] - flux_top[dim]);
                reconstruction[i][j][dim][3] += half_dt_dy * (flux_bot[dim] - flux_top[dim]);
            }
        }
    }
}

/* Calculate the fluxes in the x-direction */
Vec3D MUSCLHancockScheme::calculate_fluxes_x(const std::shared_ptr<MeshBlock>& mesh_block,
                                             const Vec4D& evolved_reconstruction)
{
    unsigned int x_cells = mesh_block->get_x_cells();
    unsigned int y_cells = mesh_block->get_y_cells();
    unsigned int dimension = mesh_block->get_dimension();
   
    Vec3D x_fluxes(x_cells + 1, Vec2D(y_cells, Vec1D(dimension, 0.0)));

    #pragma omp parallel for collapse(2)
    for (unsigned int i=0; i<x_cells + 1; i++) {
        for (unsigned int j=0; j<y_cells; j++) {
            Vec1D right_extrap(dimension, 0.0);
            Vec1D left_extrap(dimension, 0.0);

            for (unsigned int dim=0; dim<dimension; dim++) { 
                right_extrap[dim] = evolved_reconstruction[i+1][j+2][dim][1];
                left_extrap[dim] = evolved_reconstruction[i+2][j+2][dim][0];
            }

            Vec1D rs_flux = riemann->operator()(right_extrap, left_extrap, Direction::x);
        
            for (unsigned int dim=0; dim<dimension; dim++) {
                x_fluxes[i][j][dim] = rs_flux[dim];
            }
        }
    }

    return x_fluxes;
}

/* Calculate the fluxes in the y-direction */
Vec3D MUSCLHancockScheme::calculate_fluxes_y(const std::shared_ptr<MeshBlock>& mesh_block,
                                             const Vec4D& evolved_reconstruction)
{
    unsigned int x_cells = mesh_block->get_x_cells();
    unsigned int y_cells = mesh_block->get_y_cells();
    unsigned int dimension = mesh_block->get_dimension();
   
    Vec3D y_fluxes(x_cells, Vec2D(y_cells + 1, Vec1D(dimension, 0.0)));

    #pragma omp parallel for collapse(2)
    for (unsigned int i=0; i<x_cells; i++) {
        for (unsigned int j=0; j<y_cells + 1; j++) {
            Vec1D top_extrap(dimension, 0.0);
            Vec1D bottom_extrap(dimension, 0.0);

            for (unsigned int dim=0; dim<dimension; dim++) { 
                top_extrap[dim] = evolved_reconstruction[i+2][j+1][dim][3];
                bottom_extrap[dim] = evolved_reconstruction[i+2][j+2][dim][2];
            }

            Vec1D rs_flux = riemann->operator()(top_extrap, bottom_extrap, Direction::y);
            
            for (unsigned int dim=0; dim<dimension; dim++) { 
                y_fluxes[i][j][dim] = rs_flux[dim];
            }
        }
    }

    return y_fluxes;
}

/* Calculate the spatial update via the fluxes */
Vec3D MUSCLHancockScheme::calculate_spatial_update(const std::shared_ptr<MeshBlock>& mesh_block,
                                                   const Vec3D& x_fluxes,
                                                   const Vec3D& y_fluxes)
{
    double one_dx = 1.0 / mesh_block->get_dx();
    double one_dy = 1.0 / mesh_block->get_dy();

    unsigned int total_x_cells = mesh_block->get_total_x_cells();
    unsigned int total_y_cells = mesh_block->get_total_y_cells();
    unsigned int x_cells = mesh_block->get_x_cells();
    unsigned int y_cells = mesh_block->get_y_cells();
    unsigned int dimension = mesh_block->get_dimension();

    Vec3D fields_update(total_x_cells, Vec2D(total_y_cells, Vec1D(dimension, 0.0)));

    #pragma omp parallel for collapse(2)
    for (unsigned int i=0; i<x_cells; i++) {
        for (unsigned int j=0; j<y_cells; j++) {
            Vec1D cell_update = (
                (one_dx * (x_fluxes[i][j] - x_fluxes[i+1][j])) + 
                (one_dy * (y_fluxes[i][j] - y_fluxes[i][j+1]))
            );

            fields_update[i+2][j+2] = cell_update;
        }
    }

    return fields_update;
}

}

#endif
