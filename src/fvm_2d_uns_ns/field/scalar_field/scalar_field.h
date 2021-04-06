#ifndef __SCALAR_FIELD_H
#define __SCALAR_FIELD_H

#include <memory>
#include <string>

#include "../../../mesh/mesh/mesh.h"

namespace HyperFlow {

class ScalarField {

    public:

        /* Constructor */
        ScalarField();

        /* Constructor with name information */
        ScalarField(std::shared_ptr<Mesh>& _mesh,
                    const std::string& _symbol,
                    const std::string& _name,
                    const Vec1D& _initial_flow_values);

        /* Destructor */
        virtual ~ScalarField();

        /* Obtain a single flow value for a particular cell index */
        double get_flow_value(const unsigned int _cell_index);

        /* Set a single flow value for a particular cell index */
        void set_flow_value(const unsigned int _cell_index,
                            const double _flow_value);

        /* Obtain the entire collection of flow values for this field */
        Vec1D& get_full_field_flow_values();

        /* Set the entire collection of flow values for this field */
        void set_full_field_flow_values(const Vec1D& _flow_values);
        
        /* Calculate flow value gradient via least-squares estimation
         * for a particular cell index */
        double gradient(const unsigned int _cell_index);

    private:

        /* Pointer to the mesh for gradient calculations */
        std::shared_ptr<Mesh> mesh;

        /* The mathematical symbol representing the scalar field */
        std::string& symbol;

        /* The name of the scalar field */
        std::string& name;

        /* The flow value storage */
        Vec1D flow_values;

};

}

#endif
