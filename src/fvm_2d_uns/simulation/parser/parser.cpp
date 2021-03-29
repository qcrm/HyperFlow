#ifndef __PARSER_CPP
#define __PARSER_CPP

#include "parser.h"

namespace HyperFlow {

/* Constructor */
JSONParser::JSONParser()
{}

/* Constructor with JSON config */
JSONParser::JSONParser(const json& _config)
:
    config(_config)
{
    /* Start time of the simulation */
    t_start = config["simulation"]["t_start"];

    /* End time of the simulation */
    t_end = config["simulation"]["t_end"];

    /* The duration of the animation in seconds */
    anim_duration = config["simulation"]["animation"]["duration"];

    /* The frames-per-second (fps) of the animation */
    anim_fps = config["simulation"]["animation"]["fps"];

    /* Left x-coordinate of the enclosing domain */
    x_left = config["domain"]["x_left"];

    /* Right x-coordinate of the enclosing domain */
    x_right = config["domain"]["x_right"];
    
    /* Bottom y-coordinate of the enclosing domain */
    y_bottom = config["domain"]["y_bottom"];
    
    /* Top y-coordinate of the enclosing domain */
    y_top = config["domain"]["y_top"];

    /* CFL number for the explicit time-stepping */
    cfl = config["simulation"]["cfl"];
}

/* Destructor */
JSONParser::~JSONParser()
{}

/* Read access to the start time of the simulation */
const double JSONParser::get_t_start() const
{
    return t_start;
}

/* Read access to the end time of the simulation */
const double JSONParser::get_t_end() const
{
    return t_end;
}

/* Read access to the duration of the animation in seconds*/
const double JSONParser::get_anim_duration() const
{
    return anim_duration;
}

/* Read acces to the frames-per-second (fps) of the animation */
const double JSONParser::get_anim_fps() const
{
    return anim_fps;
}

/* Read access to left side x-coordinate */
const double JSONParser::get_x_left() const
{
    return x_left;
}

/* Read access to right side x-coordinate */
const double JSONParser::get_x_right() const
{
    return x_right;
}

/* Read access to bottom side y-coordinate */
const double JSONParser::get_y_bottom() const
{
    return y_bottom;
}

/* Read access to top side y-coordinate */
const double JSONParser::get_y_top() const
{
    return y_top;
}

/* Read access to the CFL number for the explicit time-stepping */
const double JSONParser::get_cfl() const
{
    return cfl;
}

/* Create model equations system */
std::shared_ptr<Model> JSONParser::generate_model()
{
    std::string model_type = config["model"]["model"];   
    if (model_type == "EulerIdealGas") {
        double gamma = config["model"]["gamma"];
        model = std::make_shared<EulerIdealGasModel>(gamma);
    } else {
        std::cout << "Model type '" << model_type << "' not supported. Returning null model." << std::endl;
    }

    return model;
}

/* Create model data output */
std::shared_ptr<DataOutput> JSONParser::generate_data_output()
{
    std::shared_ptr<EulerIdealGasDataOutput> data_output = nullptr;

    std::string model_type = config["model"]["model"];
    if (model_type == "EulerIdealGas") {
        data_output = std::make_shared<EulerIdealGasDataOutput>(
            model,
            x_left,
            x_right,
            y_bottom,
            y_top,
            t_start,
            t_end
        );
    } else {
        std::cout << "Model type '" << model_type << "' not supported. Returning null data output." << std::endl;
    }

    return data_output;
}

/* Create mesh */
std::shared_ptr<Mesh> JSONParser::generate_mesh()
{
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
   
	return mesh; 
}

/* Create initial condition */
std::shared_ptr<InitialCondition> JSONParser::generate_initial_condition()
{
    std::string initial_condition_type = config["initial_condition"]["type"];
    std::string variable_type = config["initial_condition"]["variables"];

    std::shared_ptr<InitialCondition> initcon = nullptr;

    if (initial_condition_type == "Constant") {
        initcon = generate_constant_initial_condition();
    } else if (initial_condition_type == "SphericalRiemannProblem") {
        initcon = generate_spherical_riemann_problem_initial_condition();
    } else {
        std::cout << "Initial condition '" << initial_condition_type << "' not supported. Returning null initial condition." << std::endl;
        return nullptr;
    }

    return initcon;
}

/* Create Riemann solver flux scheme */
std::shared_ptr<RiemannSolver> JSONParser::generate_riemann_solver()
{
    std::shared_ptr<RiemannSolver> riemann = nullptr;

    std::string model_type = config["model"]["model"];
    std::string riemann_type = config["scheme"]["riemann"];
    if ((model_type == "EulerIdealGas") && (riemann_type == "HLLCEuler")) {
        std::shared_ptr<EulerIdealGasModel> euler_model = std::static_pointer_cast<EulerIdealGasModel>(model);
        riemann = std::make_shared<HLLCEulerRiemannSolver>(euler_model);
    } else {
        std::cout << "Riemann solver type '" << riemann_type 
                  << "' not supported for model '" << model_type 
                  << "'. Returning null Riemann solver." << std::endl;
    }

    return riemann;
}

/* Create explicit time step calculation */
std::shared_ptr<TimeStep> JSONParser::generate_time_step()
{
    return std::make_shared<TimeStep>(model, cfl);
}

/* Create the spatial discretisation scheme */
std::shared_ptr<Scheme> JSONParser::generate_spatial_scheme(const std::shared_ptr<Model>& model,
                                                            const std::shared_ptr<RiemannSolver>& riemann,
                                                            const std::shared_ptr<TimeStep>& time_step,
                                                            const std::shared_ptr<Mesh>& mesh)
{
	std::shared_ptr<Vec1D> inlet_values = nullptr;
	Vec1D init_inlet_values;
	for (auto& elem : config["mesh"]["boundaries"]["inlet"]) {
		init_inlet_values.push_back(elem);
	}
	init_inlet_values = model->prim_to_cons(init_inlet_values);
	inlet_values = std::make_shared<Vec1D>(init_inlet_values);
    
    auto bcs = std::make_shared<BoundaryConditions>(model, inlet_values);

    std::string scheme_type = config["scheme"]["type"];
    std::shared_ptr<Scheme> scheme = nullptr;
    if (scheme_type == "Godunov") {
        scheme = std::make_shared<GodunovScheme>(model, riemann, bcs, mesh);
    } else {
        std::cout << "Scheme type '" << scheme_type << "' not supported. Exiting." << std::endl;
        return 0;
    }    
    
	return scheme;
}

/* Create ordinary differential equation time integrator */
std::shared_ptr<ODESolver> JSONParser::generate_ode_solver() {
    std::shared_ptr<ODESolver> ode_solver = nullptr;

    std::string ode_type = config["scheme"]["ode"];
    if (ode_type == "ForwardEuler") {
        ode_solver = std::make_shared<ForwardEulerODESolver>();
    } else {
        std::cout << "ODE Solver type '" << ode_type << "' not supported. Returning null ODE solver." << std::endl;
    }

    return ode_solver;
}

/* Create constant initial condition */
std::shared_ptr<ConstantInitialCondition> JSONParser::generate_constant_initial_condition()
{
    Vec1D init_state;
    for (auto& elem : config["initial_condition"]["init_state"]) {
        init_state.push_back(elem);
    }
    init_state = model->prim_to_cons(init_state);

    return std::make_shared<ConstantInitialCondition>(
        x_left, x_right, y_bottom, y_top, init_state
    );
}

/* Create spherical riemann problem initial condition */
std::shared_ptr<SphericalRiemannProblemInitialCondition> JSONParser::generate_spherical_riemann_problem_initial_condition()
{
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

    init_sphere_state.output();
    init_ext_state.output();

    init_sphere_state = model->prim_to_cons(init_sphere_state);
    init_ext_state = model->prim_to_cons(init_ext_state);

    return std::make_shared<SphericalRiemannProblemInitialCondition>(
        x_left, x_right, y_bottom, y_top, radius, x_origin, y_origin,
        init_sphere_state, init_ext_state
    );
}

}

#endif
