#ifndef __REFLECTIVE_BOUNDARY_H
#define __REFLECTIVE_BOUNDARY_H

#include "../boundary/boundary.h"

namespace HyperFlow {

class ReflectiveBoundaryCondition
:
    public BoundaryCondition
{
    
    public:
        
        /* Constructor */
        ReflectiveBoundaryCondition();

        /* Mesh block constructor */
        ReflectiveBoundaryCondition(const std::shared_ptr<MeshBlock>& _mesh_block);
        
        /* Destructor */
        virtual ~ReflectiveBoundaryCondition();

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
