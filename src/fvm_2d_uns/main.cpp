#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>

#include <nlohmann/json.hpp>

#include "scheme/ode/ode/ode.h"
#include "model/euler_ideal_gas/euler_ideal_gas.h"
#include "scheme/riemann/riemann_hllc/riemann_hllc.h"
#include "simulation/output/output.h"
#include "simulation/initcon/initcon/initcon.h"
#include "simulation/initcon/constant/constant.h"
#include "simulation/initcon/sphericalrp/sphericalrp.h"
#include "mesh/mesh/mesh.h"
#include "mesh/mesh_structured/mesh_structured.h"
#include "scheme/timestep/timestep.h"
#include "scheme/scheme/scheme.h"
#include "scheme/boundary/boundary.h"
#include "scheme/godunov/godunov.h"
#include "scheme/ode/ode/ode.h"
#include "scheme/ode/forward_euler/forward_euler.h"
#include "simulation/simulation/simulation.h"

using json = nlohmann::json;

using namespace HyperFlow;

int main(int argc, char* argv[]) {
	std::ifstream json_file(argv[1]);
    json config;
    json_file >> config;

	// Simulation
    // ==========
    double cfl = config["simulation"]["cfl"];
    double t_start = config["simulation"]["t_start"];
    double t_end = config["simulation"]["t_end"];
    double anim_duration = config["simulation"]["animation"]["duration"];
    double anim_fps = config["simulation"]["animation"]["fps"];

    // Domain 
    // ======
    double x_left = config["domain"]["x_left"];
    double x_right = config["domain"]["x_right"];
    double y_bottom = config["domain"]["y_bottom"];
    double y_top = config["domain"]["y_top"];

    // Model
    // =====
    std::string model_type = config["model"]["model"];
 
    std::shared_ptr<EulerIdealGasModel> model = nullptr;
    std::shared_ptr<EulerIdealGasDataOutput> data_output = nullptr;
    if (model_type == "EulerIdealGas") {
        model = std::make_shared<EulerIdealGasModel>(1.4);
        data_output = std::make_shared<EulerIdealGasDataOutput>(
            model,
            x_left,
            x_right,
            y_bottom,
            y_top,
            t_start,
            t_end);
    } else {
        std::cout << "Model type '" << model_type << "' not supported. Exiting." << std::endl;
        return 0;
    }

    // Mesh
    // ====
    std::shared_ptr<Mesh> mesh = nullptr;
    std::string mesh_type = config["mesh"]["type"];
    if (mesh_type == "hyperflow_structured") {
        std::string mesh_filename = config["mesh"]["file"];
        std::string base_path = "";
        std::string full_mesh_filename = base_path.append(mesh_filename);
        LoadStructuredMesh hmsh(full_mesh_filename);
        mesh = hmsh.load_structured_mesh();
    } else {
        std::cout << "Mesh type '" << mesh_type << "' not supported. Exiting." << std::endl;
        return 0;
    }

    // Initial Condition
    // =================
    std::string initial_condition_type = config["initial_condition"]["type"];
    std::string variable_type = config["initial_condition"]["variables"];
    
    std::shared_ptr<InitialCondition> initcon = nullptr;
    if (initial_condition_type == "Constant") {
        Vec1D init_state;
        
        for (auto& elem : config["initial_condition"]["init_state"]) {
            init_state.push_back(elem);
        }

        init_state = model->prim_to_cons(init_state);

        initcon = std::make_shared<ConstantInitialCondition>(
            x_left, x_right, y_bottom, y_top, t_start, t_end, init_state
        );
    }  else if (initial_condition_type == "SphericalRiemannProblem") {
        double radius = config["initial_condition"]["radius"];
        double x_origin = config["initial_condition"]["x_origin"];
        double y_origin = config["initial_condition"]["y_origin"];

        Vec1D init_sphere_state;
        Vec1D init_ext_state;
        
        for (auto& elem : config["initial_condition"]["init_sphere_state"]) {
            init_sphere_state.push_back(elem);
        }
        for (auto& elem : config["initial_condition"]["init_ext_state"]) {
            init_ext_state.push_back(elem);
        }

        init_sphere_state = model->prim_to_cons(init_sphere_state);
        init_ext_state = model->prim_to_cons(init_ext_state);

        initcon = std::make_shared<SphericalRiemannProblemInitialCondition>(
            x_left, x_right, y_bottom, y_top, radius, x_origin, y_origin,
            t_start, t_end, init_sphere_state, init_ext_state
        );
    } else {
        std::cout << "Initial condition '" << initial_condition_type << "' not supported. Exiting." << std::endl;
        return 0;
    }

    // Riemann Solver
    // ==============
    std::shared_ptr<RiemannSolver> riemann = std::make_shared<HLLCEulerRiemannSolver>(model);

    // Scheme
    // ======
    auto time_step = std::make_shared<TimeStep>(model, cfl);
	std::shared_ptr<Vec1D> inlet_values = nullptr;

	Vec1D init_inlet_values;
	for (auto& elem : config["mesh"]["boundaries"]["inlet"]) {
		init_inlet_values.push_back(elem);
	}
	init_inlet_values = model->prim_to_cons(init_inlet_values);
	inlet_values = std::make_shared<Vec1D>(init_inlet_values);
    
    auto bcs = std::make_shared<BoundaryConditions>(inlet_values);

    std::string scheme_type = config["scheme"]["type"];
    std::shared_ptr<Scheme> scheme = nullptr;
    if (scheme_type == "Godunov") {
        scheme = std::make_shared<GodunovScheme>(model, riemann, bcs, mesh);
    } else {
        std::cout << "Scheme type '" << scheme_type << "' not supported. Exiting." << std::endl;
        return 0;
    }

    // ODE Solver
    // ==========

    std::string ode_type = config["scheme"]["ode"];
    std::shared_ptr<ODESolver> ode_solver = nullptr;

    if (ode_type == "ForwardEuler") {
        ode_solver = std::make_shared<ForwardEulerODESolver>();
    } else {
        std::cout << "ODE Solver type '" << ode_type << "' not supported. Exiting." << std::endl;
        return 0;
    }

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
