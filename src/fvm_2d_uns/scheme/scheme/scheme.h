#ifndef __SCHEME_H
#define __SCHEME_H

#include "../../../share/tensor/tensor.h"
#include "../../mesh/mesh/mesh.h"

namespace HyperFlow {

class Scheme {

    public:

        /* Constructor */
        Scheme();

        /* Destructor */
        virtual ~Scheme();

        /* Apply the scheme to a mesh block */
        virtual Vec2D operator() () = 0;

};

}

#endif
