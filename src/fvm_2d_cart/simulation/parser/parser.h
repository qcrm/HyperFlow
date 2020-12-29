#ifndef __PARSER_H
#define __PARSER_H

#include <memory>
#include <vector>

#include <nlohmann/json.hpp>

#include "../../tensor/tensor.h"
#include "../../model/model/model.h"
#include "../../model/euler_ideal_gas/euler_ideal_gas.h"
#include "../../mesh/mesh_block/mesh_block.h"
#include "../../mesh/mesh/mesh.h"
#include "../../scheme/boundary/boundary/boundary.h"
#include "../../scheme/boundary/reflective/reflective_boundary.h"
#include "../../scheme/boundary/transmissive/transmissive_boundary.h"
#include "../../scheme/boundary/dirichlet/dirichlet_boundary.h"
#include "../../scheme/boundary/internal/internal_boundary.h"
#include "../initcon/initcon/initcon.h"
#include "../initcon/constant/constant.h"
#include "../initcon/doublerp/doublerp.h"
#include "../initcon/quadrp/quadrp.h"
#include "../initcon/sphericalrp/sphericalrp.h"
#include "../../scheme/riemann/riemann_hllc/riemann_hllc.h"
#include "../../scheme/timestep/timestep.h"
#include "../../scheme/godunov/godunov.h"
#include "../../scheme/ode/forward_euler/forward_euler.h"
#include "../output/euler_data_output/euler_data_output.h"

using json = nlohmann::json;

namespace HyperFlow {

class JSONParser
{
    
    public:
        
        /* Constructor */
        JSONParser();

        /* Constructor with model and JSON config */
        JSONParser(const json& _config);

        /* Destructor */
        virtual ~JSONParser();

        /* Read access to the start time of the simulation */
        const double get_t_start() const;

        /* Read access to the end time of the simulation */
        const double get_t_end() const;

        /* Read access to the duration of the animation in seconds*/
        const double get_anim_duration() const;

        /* Read acces to the frames-per-second (fps) of the animation */
        const double get_anim_fps() const;

        /* Read access to left side x-coordinate */
        const double get_x_left() const;

        /* Read access to right side x-coordinate */
        const double get_x_right() const;

        /* Read access to bottom side y-coordinate */
        const double get_y_bottom() const;

        /* Read access to top side y-coordinate */
        const double get_y_top() const;

        /* Read access to the CFL number for the explicit time-stepping */
        const double get_cfl() const;

        /* Create model equations system */
        std::shared_ptr<Model> generate_model();

        /* Create model data output */
        std::shared_ptr<DataOutput> generate_data_output();

        /* Create mesh */
        std::shared_ptr<Mesh> generate_mesh();

        /* Create initial condition */
        std::shared_ptr<InitialCondition> generate_initial_condition();

        /* Create Riemann solver flux scheme */
        std::shared_ptr<RiemannSolver> generate_riemann_solver();

        /* Create explicit time step calculation */
        std::shared_ptr<TimeStep> generate_time_step();

        /* Create the spatial discretisation scheme */
        std::shared_ptr<Scheme> generate_spatial_scheme(const std::shared_ptr<Model>& model,
                                                        const std::shared_ptr<RiemannSolver>& riemann,
                                                        const std::shared_ptr<TimeStep>& time_step);

        /* Create ordinary differential equation time integrator */
        std::shared_ptr<ODESolver> generate_ode_solver(); 

    private:

        /* Start time of the simulation */
        double t_start;

        /* End time of the simulation */
        double t_end;

        /* The duration of the animation in seconds */
        double anim_duration;

        /* The frames-per-second (fps) of the animation */
        double anim_fps;

        /* Left x-coordinate of the enclosing domain */
        double x_left;

        /* Right x-coordinate of the enclosing domain */
        double x_right;

        /* Bottom y-coordinate of the enclosing domain */
        double y_bottom;

        /* Top y-coordinate of the enclosing domain */
        double y_top;

        /* CFL number for the explicit time-stepping */
        double cfl;

        /* Pointer to the model for boundary condition
         * generation */
        std::shared_ptr<Model> model;

    	/* Config JSON */
    	json config;

        /* Create single mesh block */
        std::shared_ptr<MeshBlock> generate_mesh_block(json& sub_config);

        /* Create constant initial condition */
        std::shared_ptr<ConstantInitialCondition> generate_constant_initial_condition();

        /* Create double riemann problem initial condition */
        std::shared_ptr<DoubleRiemannProblemInitialCondition> generate_double_riemann_problem_initial_condition();

        /* Create quad riemann problem initial condition */
        std::shared_ptr<QuadRiemannProblemInitialCondition> generate_quad_riemann_problem_initial_condition();

        /* Create spherical riemann problem initial condition */
        std::shared_ptr<SphericalRiemannProblemInitialCondition> generate_spherical_riemann_problem_initial_condition();

        /* Generate an individual boundary condition from the
         * provided boundary condiguration JSON */
        std::shared_ptr<BoundaryCondition> generate_boundary_condition(
            json& boundary_config,
            std::vector<std::shared_ptr<MeshBlock> > mesh_blocks,
            std::shared_ptr<MeshBlock> mesh_block
        );

        /* Generate all boundary conditions for an individual mesh block */
        std::vector<std::shared_ptr<BoundaryCondition> > generate_block_boundary_conditions(
            json& boundaries_config,
            std::vector<std::shared_ptr<MeshBlock> > mesh_blocks,
            std::shared_ptr<MeshBlock> mesh_block
        );

};

}

#endif