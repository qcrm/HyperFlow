#ifndef __MESH_CPP
#define __MESH_CPP

#include "mesh.h"

namespace HyperFlow {

/* Constructor */
Mesh::Mesh()
{}

/* Destructor */
Mesh::~Mesh()
{}

/* Constructor with mesh information */
Mesh::Mesh(const EdgeVec1D& _edges,
           const CellVec1D& _cells,
           const std::vector<std::vector<int> >& _cell_edge_indices)
:
    edges(_edges),
    cells(_cells),
    cell_edge_indices(_cell_edge_indices)
{}

/* Get the flow values for all cells */
Vec2D Mesh::get_cells_flow_values()
{
    Vec2D cell_values = Vec2D(cells.size(), Vec1D(4, 0.0));
    for (unsigned int cell_idx=0; cell_idx<cells.size(); cell_idx++) {
        cell_values[cell_idx] = cells[cell_idx].get_flow_values();
    }
    return cell_values;
}

/* Set the flow values for all cells */
void Mesh::set_cells_flow_values(const Vec2D& cell_values)
{
    for (unsigned int cell_idx=0; cell_idx<cells.size(); cell_idx++) {
        cells[cell_idx].set_flow_values(cell_values[cell_idx]);
    }
}

/* Initialise field values */
void Mesh::initialise_field_values(const std::shared_ptr<InitialCondition>& init_con)
{
    for (unsigned int cell_idx=0; cell_idx<cells.size(); cell_idx++) {
        Vec1D centroid = cells[cell_idx].get_centroid();
        Vec1D init_flow_values = init_con->operator()(centroid[0], centroid[1]);
        cells[cell_idx].set_flow_values(init_flow_values);
    }
}

}

#endif
