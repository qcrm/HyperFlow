#ifndef __PARSER_H
#define __PARSER_H

#include <memory>
#include <vector>

#include <nlohmann/json.hpp>

#include "../../../share/tensor/tensor.h"
#include "../../model/model/model.h"
#include "../../model/euler_ideal_gas/euler_ideal_gas.h"
#include "../../mesh/mesh/mesh.h"
#include "../../mesh/mesh_structured/mesh_structured.h"
#include "../../scheme/boundary/boundary.h"
#include "../initcon/initcon/initcon.h"
#include "../initcon/constant/constant.h"
#include "../initcon/sphericalrp/sphericalrp.h"
#include "../../scheme/riemann/riemann_hllc/riemann_hllc.h"
#include "../../scheme/timestep/timestep.h"
#include "../../scheme/godunov/godunov.h"
#include "../../scheme/ode/forward_euler/forward_euler.h"
#include "../output/data_output/data_output.h"
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
                                                        const std::shared_ptr<TimeStep>& time_step,
                                                        const std::shared_ptr<Mesh>& mesh);

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

        /* Create constant initial condition */
        std::shared_ptr<ConstantInitialCondition> generate_constant_initial_condition();

        /* Create spherical riemann problem initial condition */
        std::shared_ptr<SphericalRiemannProblemInitialCondition> generate_spherical_riemann_problem_initial_condition();

};

}

#endif
