#ifndef __MUSCL_HANCOCK_H
#define __MUSCL_HANCOCK_H

#include <omp.h>

#include "../../tensor/tensor.h"
#include "../../model/model/model.h"
#include "../../scheme/riemann/riemann/riemann.h"
#include "../../mesh/mesh/mesh.h"
#include "../scheme/scheme.h"
#include "../../scheme/limiter/limiter/limiter.h"
#include "../timestep/timestep.h"

namespace HyperFlow {

class MUSCLHancockScheme
: 
    public Scheme
{

public:

    /* Constructor */
    MUSCLHancockScheme();

    /* Constructor with model equations, Riemann solver, time-step, limiter */
    MUSCLHancockScheme(const std::shared_ptr<Model> _model,
                       const std::shared_ptr<RiemannSolver> _riemann,
                       const std::shared_ptr<TimeStep> _timestep,
                       const std::shared_ptr<Limiter> _limiter,
                       const std::shared_ptr<Mesh> _mesh);

    /* Destructor */
    virtual ~MUSCLHancockScheme();

    /* Carry out the MUSCL-Hancock scheme on the provided mesh block */
    virtual Vec3D operator() (const std::shared_ptr<MeshBlock>& mesh_block);

    /* Calculate the time step via the CFL condition */
    double calculate_time_step(const std::shared_ptr<Mesh>& mesh);

    /* Number of ghost cells associated with the scheme */
    const static unsigned int ghost_cells;

private:

    /* Calculate the x-slope limiter values for each dimension */
    Vec1D calculate_slope_limiters_x(const std::shared_ptr<MeshBlock>& mesh_block,
                                     const unsigned int x_index,
                                     const unsigned int y_index);

    /* Calculate the y-slope limiter values for each dimension */
    Vec1D calculate_slope_limiters_y(const std::shared_ptr<MeshBlock>& mesh_block,
                                     const unsigned int x_index,
                                     const unsigned int y_index);

    /* Calculate the x-direction slopes for the linear reconstruction */
    Vec3D calculate_slopes_x(const std::shared_ptr<MeshBlock>& mesh_block);

    /* Calculate the y-direction slopes for the linear reconstruction */
    Vec3D calculate_slopes_y(const std::shared_ptr<MeshBlock>& mesh_block);

    /* Calculate the linear boundary extrapolations */
    Vec4D calculate_boundary_extrapolations(const std::shared_ptr<MeshBlock>& mesh_block,
                                            const Vec3D& slopes_x,
                                            const Vec3D& slopes_y);

    /* Predictor-corrector half time-step estimate of update */
    void evolve_boundary_reconstruction(const std::shared_ptr<MeshBlock>& mesh_block,
                                        Vec4D& reconstruction);

    /* Calculate the fluxes in the x-direction */
    Vec3D calculate_fluxes_x(const std::shared_ptr<MeshBlock>& mesh_block,
                             const Vec4D& evolved_reconstruction);

    /* Calculate the fluxes in the y-direction */
    Vec3D calculate_fluxes_y(const std::shared_ptr<MeshBlock>& mesh_block,
                             const Vec4D& evolved_reconstruction);

    /* Calculate the spatial update via the fluxes */
    Vec3D calculate_spatial_update(const std::shared_ptr<MeshBlock>& mesh_block,
                                   const Vec3D& x_fluxes,
                                   const Vec3D& y_fluxes);

    /* Pointer to the model equations */
    std::shared_ptr<Model> model;

    /* Pointer to the Riemann solver */
    std::shared_ptr<RiemannSolver> riemann;

    /* Pointer to the time step */
    std::shared_ptr<TimeStep> time_step;

    /* Pointer to the slope limiter */
    std::shared_ptr<Limiter> limiter;

    /* Pointer to the mesh */
    std::shared_ptr<Mesh> mesh;

    double tol;
    double omega;

};

}

#endif
