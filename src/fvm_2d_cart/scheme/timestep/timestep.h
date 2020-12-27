#ifndef __TIMESTEP_H
#define __TIMESTEP_H

#include <cmath>
#include <omp.h>

#include "../../model/model/model.h"
#include "../../mesh/mesh/mesh.h"

namespace HyperFlow {

class TimeStep {

    public:
    
        /* Constructor */
        TimeStep();
        
        /* Construct the time step with the supplied CFL
         * and pointer to the model equations */
        TimeStep(const std::shared_ptr<Model>& _model,
                 const double _cfl);
       
        /* Destructor */
        virtual ~TimeStep();
        
        /* Calculate the timestep from the provided meshblock
         * and the CFL condition */
        virtual double operator() (const std::shared_ptr<Mesh>& mesh);

    private:
       
        /* Initial CFL multiplier to stabilise
         * first few iterations */
        static double initial_cfl_multiplier;

        /* Number of steps to restrict CFL for */
        static unsigned int initial_cfl_steps;
        
        /* Pointer to the model equations */
        std::shared_ptr<Model> model;

        /* CFL value to restrict timestep */
        double cfl;

        /* Iteration counter for use in determining
         * whether to apply CFL multiplier */
        unsigned int steps;

};

}

#endif
