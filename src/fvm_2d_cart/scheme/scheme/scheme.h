#ifndef __SCHEME_H
#define __SCHEME_H

#include "../../tensor/tensor.h"
#include "../../mesh/mesh/mesh.h"

namespace HyperFlow {

class Scheme {

    public:

        /* Constructor */
        Scheme();

        /* Destructor */
        virtual ~Scheme();

        /* Apply the scheme to a mesh block */
        virtual Vec3D operator() (const std::shared_ptr<MeshBlock>& mesh_block) = 0;

        /* Calculate the time-step of the scheme */
        virtual double calculate_time_step(const std::shared_ptr<Mesh>& mesh) = 0;

};

}

#endif
