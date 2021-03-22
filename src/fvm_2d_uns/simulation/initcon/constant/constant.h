#ifndef __CONSTANT_INITCON_H
#define __CONSTANT_INITCON_H

#include "../initcon/initcon.h"

namespace HyperFlow {

class ConstantInitialCondition
:
    public InitialCondition
{
    
    public:
        
        /* Constructor */
        ConstantInitialCondition();
        
        /* Constructor with constant initial conservative flow state */
        ConstantInitialCondition(const double _x_left,
                                 const double _x_right,
                                 const double _y_bottom,
                                 const double _y_top,
                                 const double _t_start,
                                 const double _t_end,
                                 const Vec1D& _cons_init_state);
            
        /* Destructor */
        virtual ~ConstantInitialCondition();

    protected:
        
        /* Constant conservative flow values */
        Vec1D cons_init_state;

        /* Return the function for a particular dimension that
         * determines the constant flow values */
        std::function<double (const double, const double)> field_func_dim(const unsigned int dim);
};

}

#endif
