#ifndef __TRANSMISSIVE_BOUNDARY_H
#define __TRANSMISSIVE_BOUNDARY_H

#include "../boundary/boundary.h"

namespace HyperFlow {

class TransmissiveBoundaryCondition
:
    public BoundaryCondition
{
    
    public:
        
        /* Constructor */
        TransmissiveBoundaryCondition();
        
        /* Mesh block constructor */
        TransmissiveBoundaryCondition(const std::shared_ptr<MeshBlock>& _mesh_block);

        /* Destructor */
        virtual ~TransmissiveBoundaryCondition();

        /* Apply the boundary condition to the left extent */
        void apply_left();
        
        /* Apply the boundary condition to the right extent */
        void apply_right();
        
        /* Apply the boundary condition to the bottom extent */
        void apply_bottom();
        
        /* Apply the boundary condition to the top extent */
        void apply_top();

    private:

        /* Pointer to the mesh block that the boundary
         * conditions are applied to */
        std::shared_ptr<MeshBlock> mesh_block;
};

}

#endif
