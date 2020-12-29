#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

#include <nlohmann/json.hpp>

#include "scheme/ode/ode/ode.h"
#include "model/model/model.h"
#include "scheme/riemann/riemann/riemann.h"
#include "simulation/output/data_output/data_output.h"
#include "simulation/initcon/initcon/initcon.h"
#include "mesh/mesh/mesh.h"
#include "scheme/timestep/timestep.h"
#include "scheme/scheme/scheme.h"
#include "scheme/ode/ode/ode.h"
#include "simulation/simulation/simulation.h"
#include "simulation/parser/parser.h"

using json = nlohmann::json;

using namespace HyperFlow;

int main(int argc, char* argv[]) {
    // Parse the JSON configuration file
    std::ifstream json_file(argv[1]);
    json config;
    json_file >> config;

    // Create the Parser object that obtains all values
    // from the JSON config files and creates the necessary
    // object instances for the simulation
    JSONParser parser(config);

    // Obtain initial parameters
    double cfl = parser.get_cfl();
    double t_start = parser.get_t_start();
    double t_end = parser.get_t_end();
    double anim_duration = parser.get_anim_duration();
    double anim_fps = parser.get_anim_fps();

    // Equation system model
    std::shared_ptr<Model> model = parser.generate_model();

    // Model data output
    std::shared_ptr<DataOutput> data_output = parser.generate_data_output();

    // 2D Cartesian Mesh
    std::shared_ptr<Mesh> mesh = parser.generate_mesh();

    // Initial Condition
    std::shared_ptr<InitialCondition> initcon = parser.generate_initial_condition();

    // Riemann Solver
    std::shared_ptr<RiemannSolver> riemann = parser.generate_riemann_solver();

    // Time-step for explicit scheme
    std::shared_ptr<TimeStep> time_step = parser.generate_time_step();

    // Spatial discretisation scheme
    std::shared_ptr<Scheme> scheme = parser.generate_spatial_scheme(model, riemann, time_step);

    // ODE Solver
    std::shared_ptr<ODESolver> ode_solver = parser.generate_ode_solver();

    // Instantiate and solve the simulation
    // ====================================

    auto sim = std::make_shared<Simulation>(
        model,
        data_output,
        mesh,
        initcon,
        ode_solver,
        scheme,
        riemann,
        time_step,
        cfl,
        t_start,
        t_end,
        anim_duration,
        anim_fps
    );
    sim->solve();

    return 0;
}
