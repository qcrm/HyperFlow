#ifndef __MESH_CPP
#define __MESH_CPP

#include "mesh.h"

namespace HyperFlow {

/* Constructor */
Mesh::Mesh()
{}

/* Constructor from number of cells */
Mesh::Mesh(const std::vector<std::shared_ptr<MeshBlock> >& _mesh_blocks,
           const std::vector<std::vector<std::shared_ptr<BoundaryCondition> > >& _block_bcs)
:
    mesh_blocks(_mesh_blocks),
    block_bcs(_block_bcs)
{}

/* Destructor */
Mesh::~Mesh()
{}

/* Obtain write access to the mesh blocks */
std::vector<std::shared_ptr<MeshBlock> >& Mesh::get_mesh_blocks()
{
	return mesh_blocks;
}

/* Apply boundary conditions to the individual mesh blocks */
void Mesh::apply_bcs()
{
    for (unsigned int block_idx=0; block_idx<block_bcs.size(); block_idx++) {
        block_bcs[block_idx][0]->apply_bottom();
        block_bcs[block_idx][1]->apply_right();
        block_bcs[block_idx][2]->apply_top();
        block_bcs[block_idx][3]->apply_left();
    }
}

/* Initialise field values */
void Mesh::initialise_field_values(const std::shared_ptr<InitialCondition>& init_con)
{
    for (unsigned int block_idx=0; block_idx<mesh_blocks.size(); block_idx++) {
    	mesh_blocks[block_idx]->initialise_field_values(init_con);
    }
}

}

#endif