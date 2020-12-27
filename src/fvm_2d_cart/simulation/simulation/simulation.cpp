#ifndef __SIMULATION_CPP
#define __SIMULATION_CPP

#include "simulation.h"

namespace HyperFlow {

/* Constructor */
Simulation::Simulation()
{}

/* Constructor with provided pointers to all flow solver entities */
Simulation::Simulation(
    const std::shared_ptr<Model> _model,
    const std::shared_ptr<EulerIdealGasDataOutput> _data_output,
    const std::shared_ptr<Mesh> _mesh,
    const std::shared_ptr<InitialCondition> _initial_condition,
    const std::shared_ptr<ODESolver> _ode_solver,
    const std::shared_ptr<Scheme> _scheme,
    const std::shared_ptr<RiemannSolver> _riemann,
    const std::shared_ptr<TimeStep> _time_step,
    const double _cfl,
    const double _t_start,
    const double _t_end,
    const double _anim_duration,
    const double _anim_fps
)
: 
    model(_model),
    data_output(_data_output),
    mesh(_mesh),
    initial_condition(_initial_condition),
    ode_solver(_ode_solver),
    scheme(_scheme),
    riemann(_riemann),
    time_step(_time_step),
    cfl(_cfl),
    t_start(_t_start),
    t_end(_t_end),
    anim_duration(_anim_duration),
    anim_fps(_anim_fps)
{
        mesh->initialise_field_values(initial_condition);
}

/* Destructor */
Simulation::~Simulation()
{}

/* Iterate through time to solve the flow problem */
void Simulation::solve() {
    double t = t_start;
    double dt = 0.0;
    double pct_comp = 0.0;
    time_t start = time(0);

    double interval = t_end / (anim_duration * anim_fps);
    double cur_frame_out = 0.0;
    unsigned int frame_idx = 0;

    std::vector<std::shared_ptr<MeshBlock> > mesh_blocks = mesh->get_mesh_blocks();

    for (unsigned int n=1; n<Simulation::MAX_STEPS; n++) {
        double elapsed = difftime(time(0), start);
        
        pct_comp = t * 100.0 / t_end;
        
        std::cout << "Steps: "
                  << n 
                  << ", T: " 
                  << t 
                  << ", dt: " 
                  << dt 
                  << ", complete: " 
                  << pct_comp 
                  << "%, elapsed: " 
                  << elapsed 
                  << "s" 
                  << std::endl;
        
        if (t == t_end) {
            break;
        }

        mesh->apply_bcs();

        dt = scheme->calculate_time_step(mesh);
        if (t + dt > t_end) {
            dt = t_end - t;
        }

        /* Output first results data prior to any field update */
        if (cur_frame_out == 0.0) {
            data_output->to_file(mesh, frame_idx, t);
            frame_idx++;
        }

        /* Ensure the time step does not exceed the current
         * results output time */
        if ((t + dt > cur_frame_out) && (cur_frame_out > 0.0)) {
            dt = cur_frame_out - t;
        }

        /* Run one step of the solution update */

        for (unsigned int b=0; b<mesh_blocks.size(); b++) {
            auto evolve_func = std::bind(&Scheme::operator(), scheme, mesh_blocks[b]); 
            Vec3D new_fields = ode_solver->operator()(evolve_func, mesh_blocks[b]->get_fields(), dt);
            mesh_blocks[b]->set_fields(new_fields);
        }

        /* Output an interim results file when close enough to
         * the current frame output time */
        if (fabs(t - cur_frame_out) < 1e-12) {
            data_output->to_file(mesh, frame_idx, t);
            cur_frame_out += interval;
            frame_idx++;
        }

        t += dt;
    }

    data_output->to_file(mesh, frame_idx, t);
}

}

#endif
