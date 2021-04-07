#ifndef __FIELD_COLLECTION_H
#define __FIELD_COLLECTION_H

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "../../../simulation/initcon/initcon.h"
#include "../scalar_field/scalar_field.h"

namespace HyperFlow {

class FieldCollection {

    public:

        /* Constructor */
        FieldCollection();

        /* Constructor with mesh, initial condition 
         * and field name information */
        FieldCollection(std::shared_ptr<Mesh>& _mesh,
                        std::shared_ptr<InitialCondition>& _init_con,
                        const std::vector<std::pair<std::string, std::string> >& _field_symbol_names);
                        
        /* Destructor */
        virtual ~FieldCollection();

        /* Get a particular scalar field by symbol */
        ScalarField& get_scalar_field_by_symbol(const std::string& _symbol);

        /* Get a particular scalar field by name */
        ScalarField& get_scalar_field_by_name(const std::string& _name);

        /* Get particular scalar field values by symbol */
        Vec1D& get_scalar_field_values_by_symbol(const std::string& _symbol);

        /* Get particular scalar field values by name */
        Vec1D& get_scalar_field_values_by_name(const std::string& _name);

        /* Get an individual scalar field value by symbol */
        double get_scalar_field_flow_value_by_symbol(const std::string& _symbol);

        /* Get an particular scalar field value by name */
        double get_scalar_field_flow_value_by_name(const std::string& _name);

    private:

        /* Pointer to the mesh */
        std::shared_ptr<Mesh> mesh;

        /* Pointer to the initial condition isntance */
        std::shared_ptr<InitialCondition> init_con;
       
        /* List of the field names/symbols */
        std::vector<std::pair<std::string, std::string> > field_symbol_names;

        /* Symbol to name storage mapping */
        std::unordered_map<std::string, std::string> symbol_names;

        /* Symbol to field storage mapping */
        std::unordered_map<std::string, ScalarField> fields;

        /* Populate the field name -> symbol map */
        void populate_name_to_symbol_map();

        /* Populate the map of field symbols to names
         * by creating each ScalarField */
        void populate_scalar_field_map();

};

}

#endif
