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
    std::shared_ptr<DataOutput> data_output = nullptr;

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
    /* Make mesh blocks */
    std::vector<std::shared_ptr<MeshBlock> > mesh_blocks;
    for (unsigned int block_idx=0; block_idx<config["mesh"]["blocks"].size(); block_idx++) {
        std::shared_ptr<MeshBlock> mb = generate_mesh_block(
            config["mesh"]["blocks"][block_idx]
        );
        mesh_blocks.push_back(mb);
    }

    /* Make boundary conditions */
    std::vector<std::vector<std::shared_ptr<BoundaryCondition> > > block_bcs;
    for (unsigned int block_idx=0; block_idx<config["mesh"]["blocks"].size(); block_idx++) {
        std::vector<std::shared_ptr<BoundaryCondition> > bcs = generate_block_boundary_conditions(
            config["mesh"]["blocks"][block_idx]["boundaries"], mesh_blocks, mesh_blocks[block_idx]
        );
        block_bcs.push_back(bcs);
    }

    std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>(mesh_blocks, block_bcs);
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
    } else if (initial_condition_type == "DoubleRiemannProblem") {
        initcon = generate_double_riemann_problem_initial_condition();
    } else if (initial_condition_type == "QuadRiemannProblem") {
        initcon = generate_quad_riemann_problem_initial_condition();
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
    std::shared_ptr<Scheme> scheme = nullptr;

    std::string scheme_type = config["scheme"]["type"];
    if (scheme_type == "Godunov") {
        scheme = std::make_shared<GodunovScheme>(model, riemann, time_step);
    } else if (scheme_type == "MUSCLHancock") {
        std::string limiter_type = config["scheme"]["limiter"];

        std::shared_ptr<Limiter> limiter = nullptr;

        if (limiter_type == "MinBee") {
            limiter = std::make_shared<MinBeeLimiter>(0.0, 1.0, 1.0, 1e-6);
        } else if (limiter_type == "VanLeer") {
            limiter = std::make_shared<VanLeerLimiter>(0.0, 1.0, 1.0, 1e-6);
        } else {
            std::cout << "Limiter type '" << limiter_type << "' not supported. Exiting." << std::endl;
            return 0;
        }

        scheme = std::make_shared<MUSCLHancockScheme>(model, riemann, time_step, limiter, mesh);
    } else {
        std::cout << "Scheme type '" << scheme_type << "' not supported. Returning null scheme." << std::endl;
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

/* Create single mesh block */
std::shared_ptr<MeshBlock> JSONParser::generate_mesh_block(json& sub_config)
{
    Vec1D extents {
        sub_config["extents"]["x_left"],
        sub_config["extents"]["x_right"],
        sub_config["extents"]["y_bottom"],
        sub_config["extents"]["y_top"]
    };

    unsigned int x_cells = sub_config["x_cells"];
    unsigned int y_cells = sub_config["y_cells"];

    // Set the number of ghost cells depending upon the scheme
    unsigned int ghost_cells = 1;
    std::string scheme_type = config["scheme"]["type"];
    if (scheme_type == "Godunov") {
        ghost_cells = 1;
    } else if (scheme_type == "MUSCLHancock") {
        ghost_cells = 2;
    } else {
        std::cout << "Scheme type '"
                  << scheme_type 
                  << "' not supported. Number of ghost cells (" 
                  << ghost_cells 
                  << ") may be incorrect." 
                  << std::endl;
    }

    // Set the dimension of the mesh depending upon the equation system
    unsigned int dimension = 4;
    std::string model_type = config["model"]["model"];
    if (model_type == "EulerIdealGas") {
        dimension = 4;
    } else {
        std::cout << "Model type '" 
                  << model_type
                  << "' not supported. Number of mesh dimensions ("
                  << dimension 
                  << ")may be incorrect."
                  << std::endl;
    }
    
    std::shared_ptr<MeshBlock> mb = std::make_shared<MeshBlock>(
        extents, x_cells, y_cells, ghost_cells, dimension
    );
    return mb;
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

/* Create double riemann problem initial condition */
std::shared_ptr<DoubleRiemannProblemInitialCondition> JSONParser::generate_double_riemann_problem_initial_condition()
{
    std::string variable_type = config["initial_condition"]["variables"];

    double x_lm_interface = config["initial_condition"]["x_lm_interface"];
    double x_mr_interface = config["initial_condition"]["x_mr_interface"];

    Vec1D init_left_state;
    Vec1D init_middle_state;
    Vec1D init_right_state;
    
    for (auto& elem : config["initial_condition"]["init_left_state"]) {
        init_left_state.push_back(elem);
    }
    for (auto& elem : config["initial_condition"]["init_middle_state"]) {
        init_middle_state.push_back(elem);
    }
    for (auto& elem : config["initial_condition"]["init_right_state"]) {
        init_right_state.push_back(elem);
    }

    if (variable_type == "Primitive") {
        init_left_state = model->prim_to_cons(init_left_state);
        init_middle_state = model->prim_to_cons(init_middle_state);
        init_right_state = model->prim_to_cons(init_right_state);
    }

    return std::make_shared<DoubleRiemannProblemInitialCondition>(
        x_left, x_right, y_bottom, y_top, x_lm_interface, x_mr_interface, 
        init_left_state, init_middle_state, init_right_state
    );
}

/* Create quad riemann problem initial condition */
std::shared_ptr<QuadRiemannProblemInitialCondition> JSONParser::generate_quad_riemann_problem_initial_condition()
{
    std::string variable_type = config["initial_condition"]["variables"];

    double x_interface = config["initial_condition"]["x_interface"];
    double y_interface = config["initial_condition"]["y_interface"];

    Vec1D xy_interface(2, 0.0);
    xy_interface[0] = x_interface;
    xy_interface[1] = y_interface;

    Vec1D init_nw_state;
    Vec1D init_ne_state;
    Vec1D init_sw_state;
    Vec1D init_se_state;
    
    for (auto& elem : config["initial_condition"]["init_nw_state"]) {
        init_nw_state.push_back(elem);
    }
    for (auto& elem : config["initial_condition"]["init_ne_state"]) {
        init_ne_state.push_back(elem);
    }
    for (auto& elem : config["initial_condition"]["init_sw_state"]) {
        init_sw_state.push_back(elem);
    }
    for (auto& elem : config["initial_condition"]["init_se_state"]) {
        init_se_state.push_back(elem);
    }

    if (variable_type == "Primitive") {
        init_nw_state = model->prim_to_cons(init_nw_state);
        init_ne_state = model->prim_to_cons(init_ne_state);
        init_sw_state = model->prim_to_cons(init_sw_state);
        init_se_state = model->prim_to_cons(init_se_state);
    }

    return std::make_shared<QuadRiemannProblemInitialCondition>(
        x_left, x_right, y_bottom, y_top, xy_interface, 
        init_nw_state, init_ne_state, init_sw_state, init_se_state
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

/* Generate an individual boundary condition from the
 * provided boundary condiguration JSON */
std::shared_ptr<BoundaryCondition> JSONParser::generate_boundary_condition(
    json& boundary_config,
    std::vector<std::shared_ptr<MeshBlock> > mesh_blocks,
    std::shared_ptr<MeshBlock> mesh_block
)
{
    std::shared_ptr<BoundaryCondition> bc = nullptr;
    std::string boundary_type = boundary_config["type"];

    if (boundary_type == "Dirichlet") {
        std::string variables_type = boundary_config["variables"];
        Vec1D flow_state;
        
        for (auto& elem : boundary_config["flow_state"]) {
            flow_state.push_back(elem);
        }
        
        if (variables_type == "Primitive") {
            flow_state = model->prim_to_cons(flow_state);
        }

        bc = std::make_shared<DirichletBoundaryCondition>(mesh_block, flow_state);
    } else if (boundary_type == "Reflective") {
        bc = std::make_shared<ReflectiveBoundaryCondition>(mesh_block);
    } else if (boundary_type == "Transmissive") {
        bc = std::make_shared<TransmissiveBoundaryCondition>(mesh_block);
    } else if (boundary_type == "Internal") {
        unsigned int block_idx = boundary_config["block"];
        bc = std::make_shared<InternalBoundaryCondition>(mesh_block, mesh_blocks[block_idx]);
    } else {
        std::cout << "Boundary condition '" << boundary_type << "' not recognised." << std::endl;
    }

    return bc;
}

/* Generate all boundary conditions for an individual mesh block */
std::vector<std::shared_ptr<BoundaryCondition> > JSONParser::generate_block_boundary_conditions(
    json& boundaries_config,
    std::vector<std::shared_ptr<MeshBlock> > mesh_blocks,
    std::shared_ptr<MeshBlock> mesh_block
)
{
    std::vector<std::shared_ptr<BoundaryCondition> > bcs;

    std::shared_ptr<BoundaryCondition> bc_bottom = generate_boundary_condition(boundaries_config["bottom"], mesh_blocks, mesh_block);
    std::shared_ptr<BoundaryCondition> bc_right = generate_boundary_condition(boundaries_config["right"], mesh_blocks, mesh_block);
    std::shared_ptr<BoundaryCondition> bc_top = generate_boundary_condition(boundaries_config["top"], mesh_blocks, mesh_block);
    std::shared_ptr<BoundaryCondition> bc_left = generate_boundary_condition(boundaries_config["left"], mesh_blocks, mesh_block);

    bcs.push_back(bc_bottom);
    bcs.push_back(bc_right);
    bcs.push_back(bc_top);
    bcs.push_back(bc_left);

    return bcs;
}

}

#endif