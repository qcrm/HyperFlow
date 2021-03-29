#ifndef __SPHERICAL_RIEMANN_PROBLEM_H
#define __SPHERICAL_RIEMANN_PROBLEM_H

#include "../initcon/initcon.h"

namespace HyperFlow {

class SphericalRiemannProblemInitialCondition
:
    public InitialCondition
{
    
    public:
        
        /* Constructor */
        SphericalRiemannProblemInitialCondition();
        
        /* Constructor with spatio-temporal extent, sphere
         * geometry and internal and external flow states */
        SphericalRiemannProblemInitialCondition(const double _x_left,
                                                const double _x_right,
                                                const double _y_bottom,
                                                const double _y_top,
                                                const double _radius,
                                                const double _x_origin,
                                                const double _y_origin,
                                                const Vec1D& _cons_sphere_state,
                                                const Vec1D& _cons_ext_state);
            
        /* Destructor */
        virtual ~SphericalRiemannProblemInitialCondition();

        /* Obtain the radius of the sphere flow state */
        double get_radius();
        
        /* Obtain the x coordinate of the sphere flow state centre */
        double get_x_origin();
        
        /* Obtain the y coordinate of the sphere flow state centre */
        double get_y_origin();

    protected:

        /* Radius of the sphere flow state */
        double radius;
        
        /* x coordinate of the sphere flow state centre */
        double x_origin;
        
        /* y coordinate of the sphere flow state centre */
        double y_origin;
        
        /* Conservative flow values inside the sphere */
        Vec1D cons_sphere_state;
        
        /* Conservative flow values external to the sphere */
        Vec1D cons_ext_state;

        /* Return the function for a particular dimension that
         * determines the flow values internal and external
         * to the spherical region */
        std::function<double (const double, const double)> field_func_dim(const unsigned int dim);
};

}

#endif
