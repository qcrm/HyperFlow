#ifndef __SIMULATION_H
#define __SIMULATION_H

#include <iostream>
#include <time.h>

#include "../output/data_output/data_output.h"
#include "../initcon/initcon/initcon.h"
#include "../../model/model/model.h"
#include "../../scheme/riemann/riemann/riemann.h"
#include "../../mesh/mesh/mesh.h"
#include "../../scheme/godunov/godunov.h"
#include "../../scheme/ode/ode/ode.h"
#include "../../scheme/timestep/timestep.h"

namespace HyperFlow {

class Simulation
{
    
    public:
       
        /* Constructor */
        Simulation();
        
        /* Constructor with provided pointers to all flow solver entities */
        Simulation(const std::shared_ptr<Model> _model,
                   const std::shared_ptr<DataOutput> _data_output,
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
                   const double _anim_fps);
        
        /* Destructor */
        virtual ~Simulation();
        
        /* Iterate through time to solve the flow problem */
        void solve();
        
    private:

        /* Maximum number of timestep iterations */
        static const unsigned int MAX_STEPS = 1000000000;

        /* Pointer to the model equations */
        std::shared_ptr<Model> model;
        
        /* Pointer to the data output */
        std::shared_ptr<DataOutput> data_output;
        
        /* Pointer to the mesh */
        std::shared_ptr<Mesh> mesh;
       
        /* Pointer to the initial condition */
        std::shared_ptr<InitialCondition> initial_condition;

        /* Pointer to the ODE solver */
        std::shared_ptr<ODESolver> ode_solver;
        
        /* Pointer to the spatial discretisation scheme */
        std::shared_ptr<Scheme> scheme;
        
        /* Pointer to the Riemann problem used in the spatial
         * discretisation scheme */
        std::shared_ptr<RiemannSolver> riemann;
        
        /* Pointer to the timestep calculator */
        std::shared_ptr<TimeStep> time_step;
        
        /* Courant number for CFL condition */
        double cfl;
        
        /* Start time in seconds */
        double t_start;
        
        /* End time in seconds */
        double t_end;
        
        /* Animation duration in seconds */
        double anim_duration;
        
        /* Animation frames per second */
        double anim_fps;

};

}

#endif
