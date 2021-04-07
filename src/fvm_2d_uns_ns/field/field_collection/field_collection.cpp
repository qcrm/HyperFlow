#ifndef __SCALAR_FIELD_CPP
#define __SCALAR_FIELD_CPP

#include "field_collection.h"

namespace HyperFlow {

/* Constructor */
FieldCollection::FieldCollection()
{}

/* Constructor with mesh, initial condition 
 * and field name information */
FieldCollection::FieldCollection(std::shared_ptr<Mesh>& _mesh,
                                 std::shared_ptr<InitialCondition>& _init_con,
                                 std::vector<std::pair<std::string, std::string> >& _field_symbol_names)
:
    mesh(_mesh),
    init_con(_init_con),
    field_symbol_names(_field_symbol_names)
{
    populate_name_to_symbol_map();
    populate_scalar_field_map();
}
                
/* Destructor */
FieldCollection::~FieldCollection()
{}

/* Get a particular scalar field by symbol */
ScalarField& FieldCollection::get_scalar_field_by_symbol(const std::string& _symbol)
{
   return fields[symbol]; 
}

/* Get a particular scalar field by name */
ScalarField& FieldCollection::get_scalar_field_by_name(const std::string& _name)
{
   return fields[symbol_map[_name]]; 

}

/* Get particular scalar field values by symbol */
Vec1D& FieldCollection::get_scalar_field_values_by_symbol(const std::string& _symbol)
{
    return fields[symbol].get_full_field_flow_values();
}

/* Get particular scalar field values by name */
Vec1D& FieldCollection::get_scalar_field_values_by_name(const std::string& _name)
{
    return fields[symbol_map[_name]].get_full_field_flow_values();
}

/* Get an individual scalar field value by symbol */
double FieldCollection::get_scalar_field_flow_value_by_symbol(const std::string& _symbol,
                                                              const unsigned int _cell_index)
{
    return fields[symbol].get_flow_value(_cell_index);
}

/* Get an particular scalar field value by name */
double FieldCollection::get_scalar_field_flow_value_by_name(const std::string& _name,
                                                            const unsigned int _cell_index)
{
    return fields[symbol_map[_name]].get_flow_value(_cell_index);
}

/* Populate the field name -> symbol map */
void FieldCollection::populate_name_to_symbol_map()
{
    for (const auto& symname : field_symbol_names) {
        symbol_names[symname.second] = symname.first;
    }
}

/* Populate the map of field symbols to names
 * by creating each ScalarField */
void FieldCollection::populate_scalar_field_map()
{
    for (const auto& symname : field_symbol_names) {
        // TODO: Create initial flow values
        
        fields[symname.first] = ScalarField(mesh, 
                                            symname.first,
                                            symname.second,
                                            initial_flow_values);

    }
}

}

#endif
