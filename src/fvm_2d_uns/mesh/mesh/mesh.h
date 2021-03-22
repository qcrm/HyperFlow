#ifndef __MESH_H
#define __MESH_H

#include "../../../share/tensor/tensor.h"
#include "../../geometry/edge/edge.h"
#include "../../geometry/polygon/polygon.h"
#include "../../mesh/cell/cell.h"
#include "../../simulation/initcon/initcon/initcon.h"

namespace HyperFlow {

class Mesh {

    public:

        /* Constructor */
        Mesh();

        /* Constructor with mesh information */
        Mesh(const EdgeVec1D& _edges,
             const CellVec1D& _cells,
             const std::vector<std::vector<int> >& _cell_edge_indices);

        /* Destructor */
        virtual ~Mesh();

        /* Initialise field values */
        void initialise_field_values(const std::shared_ptr<InitialCondition>& init_con);

        /* Get the flow values for all cells */
        Vec2D get_cells_flow_values();

        /* Set the flow values for all cells */
        void set_cells_flow_values(const Vec2D& cell_values);

    private:

        /* The vector of all edges */
        EdgeVec1D edges;

        /* The vector of all cells */
        CellVec1D cells;

        /* The vector listing all edge indices of a cell into the 'edges' EdgeVec1D */
        std::vector<std::vector<int> > cell_edge_indices;

};

}

#endif
