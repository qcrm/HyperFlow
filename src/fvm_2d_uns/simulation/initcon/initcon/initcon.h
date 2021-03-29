#ifndef __INITCON_H
#define __INITCON_H

#include <functional>
#include <memory>

#include "../../../../share/tensor/tensor.h"

namespace HyperFlow {

class InitialCondition
{

    public:
       
        /* Constructor */
        InitialCondition();
        
        /* Constructor from spatial and temporal extents */
        InitialCondition(const double _x_left,
                         const double _x_right,
                         const double _y_bottom,
                         const double _y_top,
                         const unsigned int _dimension);
            
        /* Destructor */
        virtual ~InitialCondition();

        /* Obtain the initial flow values at the
         * provided coordinates */
        virtual Vec1D operator() (const double x, const double y);
        
        /* Obtain the left extent of the simulation */
        virtual double get_x_left();
        
        /* Obtain the right extent of the simultion */
        virtual double get_x_right();
        
        /* Obtain the bottom extent of the simulation */
        virtual double get_y_bottom();
        
        /* Obtain the top extent of the simulation */
        virtual double get_y_top();
        
        /* Obtain the dimensionality of the model system */
        virtual unsigned int get_dimension();

    protected:
        
        /* Left extent of the simulation */
        double x_left;
        
        /* Right extent of the simulation */
        double x_right;
        
        /* Bottom extent of the simulation */
        double y_bottom;
        
        /* Top extent of the simulation */
        double y_top;
        
        /* Dimensionality of the model system */
        unsigned int dimension;

        /* Return the flow value function at a particular dimension */
        virtual std::function<double (const double, const double)> field_func_dim(const unsigned int dim) = 0;
};

}

# endif
