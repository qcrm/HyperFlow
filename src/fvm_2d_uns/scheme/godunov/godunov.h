#ifndef __GODUNOV_H
#define __GODUNOV_H

#include <cmath>
#include <iostream>

#include "../../../share/tensor/tensor.h"
#include "../../model/model/model.h"
#include "../../scheme/riemann/riemann/riemann.h"
#include "../../mesh/mesh/mesh.h"
#include "../scheme/scheme.h"

namespace HyperFlow {

class GodunovScheme
: 
    public Scheme
{

public:

    /* Constructor */
    GodunovScheme();

    /* Parameterised constructor */
    GodunovScheme(
        const std::shared_ptr<Model> _model,
        const std::shared_ptr<RiemannSolver> _riemann,
        std::shared_ptr<Mesh> _mesh
    );

    /* Destructor */
    virtual ~GodunovScheme();
    
    /* Carry out the Godunov scheme on the provided mesh */
    virtual Vec2D operator()();
        
private:
    std::shared_ptr<Model> model;
    std::shared_ptr<RiemannSolver> riemann;
    std::shared_ptr<Mesh> mesh;

};

}

#endif
