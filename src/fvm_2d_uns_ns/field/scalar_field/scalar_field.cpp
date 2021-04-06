#ifndef __SCALAR_FIELD_CPP
#define __SCALAR_FIELD_CPP

#include "scalar_field.h"

namespace HyperFlow {

/* Constructor */
ScalarField::ScalarField();

/* Constructor with name information */
ScalarField::ScalarField(std::shared_ptr<Mesh>& _mesh,
                         const std::string& _symbol,
                         const std::string& _name,
                         const Vec1D& _initial_flow_values)
:
    mesh(_mesh),
    symbol(_symbol),
    name(_name),
    flow_values(_initial_flow_values)
{}

/* Destructor */
ScalarField::~ScalarField()
{}

/* Obtain a single flow value for a particular cell index */
double ScalarField::get_flow_value(const unsigned int _cell_index)
{
    return flow_values[_cell_index];
}

/* Set a single flow value for a particular cell index */
void ScalarField::set_flow_value(const unsigned int _cell_index,
                                 const double _flow_value)
{
    flow_values[_cell_index] = _flow_value;
}

/* Obtain the entire collection of flow values for this field */
Vec1D& ScalarField::get_full_field_flow_values()
{
    return flow_values;
}

/* Set the entire collection of flow values for this field */
void ScalarField::set_full_field_flow_values(const Vec1D& _flow_values)
{
    flow_values = _flow_values;
}

/* Calculate flow value gradient via least-squares estimation
 * for a particular cell index */
double ScalarField::gradient(const unsigned int _cell_index)
{
    // TODO
    return 0.0;
}

};

}

#endif
