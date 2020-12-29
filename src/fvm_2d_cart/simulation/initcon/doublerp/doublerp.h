#ifndef __DOUBLE_RIEMANN_PROBLEM_H
#define __DOUBLE_RIEMANN_PROBLEM_H

#include "../initcon/initcon.h"

namespace HyperFlow {

class DoubleRiemannProblemInitialCondition
:
    public InitialCondition
{
    
    public:
        
        /* Constructor */
        DoubleRiemannProblemInitialCondition();
        
        /* Parameterised constructor */
        DoubleRiemannProblemInitialCondition(const double _x_left,
                                             const double _x_right,
                                             const double _y_bottom,
                                             const double _y_top,
                                             const double _x_lm_interface,
                                             const double _x_mr_interface,
                                             const Vec1D& _cons_left_state,
                                             const Vec1D& _cons_middle_state,
                                             const Vec1D& _cons_right_state);
            
        /* Destructor */
        virtual ~DoubleRiemannProblemInitialCondition();

    protected:
        
        /* Left-Middle interface x-coordinate */
        double x_lm_interface;

        /* Middle-Right interface x-coordinate */
        double x_mr_interface;

        /* Conservative flow values for the left state */
        Vec1D cons_left_state;

        /* Conservative flow values for the middle state */
        Vec1D cons_middle_state;

        /* Conservative flow values for the right state */
        Vec1D cons_right_state;

        /* Return the function for a particular dimension that
         * determines the flow values for each state */
        std::function<double (const double, const double)> field_func_dim(const unsigned int dim);
};

}

#endif
