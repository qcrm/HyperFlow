#ifndef __MESH_H
#define __MESH_H

#include "../../tensor/tensor.h"
#include "../../simulation/initcon/initcon/initcon.h"
#include "../../scheme/boundary/boundary/boundary.h"
#include "../mesh_block/mesh_block.h"

namespace HyperFlow {

class Mesh {

    public:

        /* Constructor */
        Mesh();

        /* Constructor from number of cells */
        Mesh(const std::vector<std::shared_ptr<MeshBlock> >& _mesh_blocks,
        	 const std::vector<std::vector<std::shared_ptr<BoundaryCondition> > >& _block_bcs);

        /* Destructor */
        virtual ~Mesh();

        /* Obtain write access to the mesh blocks */
        std::vector<std::shared_ptr<MeshBlock> >& get_mesh_blocks();

        /* Apply boundary conditions to the individual mesh blocks */
        void apply_bcs();

        /* Initialise field values */
        void initialise_field_values(const std::shared_ptr<InitialCondition>& init_con);

    private:

    	/* Collection of the mesh blocks making up the overall mesh */
    	std::vector<std::shared_ptr<MeshBlock> > mesh_blocks;
        
    	/* Two-dimensional vector of boundary conditions. Outer dimension 
    	 * is the size of number of meshes. Inner dimension represents 
    	 * bottom, right, top and left boundaries for each block. */
		std::vector<std::vector<std::shared_ptr<BoundaryCondition> > > block_bcs;

};

}

#endif