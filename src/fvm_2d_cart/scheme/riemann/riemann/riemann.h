#ifndef __RIEMANN_H
#define __RIEMANN_H

#include "../../../tensor/tensor.h" 

namespace HyperFlow {

/* Allow the Riemann solver to use x or y direction
 * for flux calculation */
enum class Direction { x, y };

class RiemannSolver
{

public:
    
    /* Constructor */
    RiemannSolver();

    /* Destructor */
    virtual ~RiemannSolver();
    
    /* Evaluate Riemann problem */
    virtual Vec1D operator() (const Vec1D& cons_left,
                              const Vec1D& cons_right,
                              const Direction& dir) = 0;
};

}

#endif
