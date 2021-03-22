#ifndef __MESH_STRUCTURED_H
#define __MESH_STRUCTURED_H

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "../../../share/tensor/tensor.h"
#include "../../geometry/edge/edge.h"
#include "../../geometry/polygon/polygon.h"
#include "../../mesh/cell/cell.h"
#include "../../mesh/mesh/mesh.h"

namespace HyperFlow {

class LoadStructuredMesh {

    public:

        /* Constructor */
        LoadStructuredMesh();

        /* Constructor with filename */
        LoadStructuredMesh(const std::string& _filename);

        /* Destructor */
        virtual ~LoadStructuredMesh();

        /* Load the HyperFlor native format structured mesh from disk */
        std::shared_ptr<Mesh> load_structured_mesh();

    private:

        /* The filename of the structured mesh */
        std::string filename;

};

}

#endif
