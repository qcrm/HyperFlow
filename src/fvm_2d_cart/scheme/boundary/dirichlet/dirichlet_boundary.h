#ifndef __DIRICHLET_BOUNDARY_H
#define __DIRICHLET_BOUNDARY_H

#include "../boundary/boundary.h"

namespace HyperFlow {

class DirichletBoundaryCondition
:
    public BoundaryCondition
{
    
    public:
        
        /* Constructor */
        DirichletBoundaryCondition();

        /* Constructor with specified conservative flow values */
        DirichletBoundaryCondition(const std::shared_ptr<MeshBlock>& _mesh_block,
                                   const Vec1D& _cons_flow_values);
        
        /* Destructor */
        virtual ~DirichletBoundaryCondition();

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

        /* Specified conservative flow values */
        Vec1D cons_flow_values;
};

}

#endif
