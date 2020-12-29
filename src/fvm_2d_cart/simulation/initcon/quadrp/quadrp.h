#ifndef __QUAD_RIEMANN_PROBLEM_H
#define __QUAD_RIEMANN_PROBLEM_H

#include "../initcon/initcon.h"

namespace HyperFlow {

class QuadRiemannProblemInitialCondition
:
    public InitialCondition
{
    
    public:
        
        /* Constructor */
        QuadRiemannProblemInitialCondition();
        
        /* Constructor with spatio-temporal extent, interface
         * location and all four flow states */
        QuadRiemannProblemInitialCondition(const double _x_left,
                                           const double _x_right,
                                           const double _y_bottom,
                                           const double _y_top,
                                           const Vec1D& _xy_interface,
                                           const Vec1D& _cons_nw_state,
                                           const Vec1D& _cons_ne_state,
                                           const Vec1D& _cons_sw_state,
                                           const Vec1D& _cons_se_state);
            
        /* Destructor */
        virtual ~QuadRiemannProblemInitialCondition();
        
        /* Obtain the x-y coordinates of the quadrant interface */
        Vec1D get_xy_interface();

    protected:

        /* x-y coordinates of the quadrant interface */
        Vec1D xy_interface;
        
        /* Conservative flow values for the north-west state */
        Vec1D cons_nw_state;

        /* Conservative flow values for the north-east state */
        Vec1D cons_ne_state;

        /* Conservative flow values for the south-west state */
        Vec1D cons_sw_state;

        /* Conservative flow values for the south-east state */
        Vec1D cons_se_state;

        /* Return the function for a particular dimension that
         * determines the flow values for each quadrant */
        std::function<double (const double, const double)> field_func_dim(const unsigned int dim);
};

}

#endif
