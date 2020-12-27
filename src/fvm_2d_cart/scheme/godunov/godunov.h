#ifndef __GODUNOV_H
#define __GODUNOV_H

#include "../../tensor/tensor.h"
#include "../../model/model/model.h"
#include "../../scheme/riemann/riemann/riemann.h"
#include "../../mesh/mesh/mesh.h"
#include "../scheme/scheme.h"
#include "../timestep/timestep.h"

namespace HyperFlow {

class GodunovScheme
: 
    public Scheme
{

public:

    /* Constructor */
    GodunovScheme();

    /* Constructor with model equations, Riemann solver, time-step */
    GodunovScheme(const std::shared_ptr<Model>& _model,
                  const std::shared_ptr<RiemannSolver>& _riemann,
                  const std::shared_ptr<TimeStep>& _timestep);

    /* Destructor */
    virtual ~GodunovScheme();

    /* Carry out the Godunov scheme on the provided mesh */
    virtual Vec3D operator() (const std::shared_ptr<MeshBlock>& mesh_block);

    /* Calculate the time step via the CFL condition */
    double calculate_time_step(const std::shared_ptr<Mesh>& mesh);

    /* Number of ghost cells associated with the scheme */
    const static unsigned int ghost_cells;

private:

    /* Calculate the fluxes in the x-direction */
    Vec3D calculate_fluxes_x(const std::shared_ptr<MeshBlock>& mesh_block);

    /* Calculate the fluxes in the y-direction */
    Vec3D calculate_fluxes_y(const std::shared_ptr<MeshBlock>& mesh_block);

    /* Calculate the spatial update via the fluxes */
    Vec3D calculate_spatial_update(const std::shared_ptr<MeshBlock>& mesh_block,
                                   const Vec3D& x_fluxes,
                                   const Vec3D& y_fluxes);

    /* Pointer to the model equations */
    std::shared_ptr<Model> model;

    /* Pointer to the Riemann solver */
    std::shared_ptr<RiemannSolver> riemann;

    /* Pointer to the time step */
    std::shared_ptr<TimeStep> timestep;

};

}

#endif
