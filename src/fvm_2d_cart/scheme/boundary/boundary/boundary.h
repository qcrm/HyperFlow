#ifndef __BOUNDARY_H
#define __BOUNDARY_H

#include <omp.h>

#include "../../../mesh/mesh_block/mesh_block.h"

namespace HyperFlow {

class BoundaryCondition {

    public:
       
        /* Constructor */
        BoundaryCondition();
       
        /* Destructor */
        virtual ~BoundaryCondition();

        /* Apply the boundary condition to the left extent */
        virtual void apply_left() = 0;
        
        /* Apply the boundary condition to the right extent */
        virtual void apply_right() = 0;
        
        /* Apply the boundary condition to the bottom extent */
        virtual void apply_bottom() = 0;
        
        /* Apply the boundary condition to the top extent */
        virtual void apply_top() = 0;
};

}

#endif
