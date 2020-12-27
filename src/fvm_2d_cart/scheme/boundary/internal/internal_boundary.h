#ifndef __INTERNAL_BOUNDARY_H
#define __INTERNAL_BOUNDARY_H

#include "../boundary/boundary.h"

namespace HyperFlow {

class InternalBoundaryCondition
:
    public BoundaryCondition
{
    
    public:
        
        /* Constructor */
        InternalBoundaryCondition();

        /* Mesh block constructor */
        InternalBoundaryCondition(const std::shared_ptr<MeshBlock>& _internal_block,
                                  const std::shared_ptr<MeshBlock>& _external_block);
        
        /* Destructor */
        virtual ~InternalBoundaryCondition();

        /* Apply the boundary condition to the left extent */
        void apply_left();
        
        /* Apply the boundary condition to the right extent */
        void apply_right();
        
        /* Apply the boundary condition to the bottom extent */
        void apply_bottom();
        
        /* Apply the boundary condition to the top extent */
        void apply_top();

    private:

        /* Pointer to the internal block that the boundary
         * conditions are applied to */
        std::shared_ptr<MeshBlock> internal_block;
        
        /* Pointer to the neighbouring block that the ghost
         * cells are copied from */
        std::shared_ptr<MeshBlock> external_block;

};

}

#endif
