#ifndef __CELL_CENTRED_FIELD_MESH_CPP
#define __CELL_CENTRED_FIELD_MESH_CPP

#include "cell_mesh.h"

namespace HyperFlow {

/* Constructor */
CellCentredFieldMesh::CellCentredFieldMesh()
{}

/* Constructor from number of cells */
CellCentredFieldMesh::CellCentredFieldMesh(const double _x_left,
                                           const double _x_right,
                                           const unsigned int _cells,
                                           const unsigned int _ghost_cells)
:
    x_left(_x_left),
    x_right(_x_right),
    cells(_cells),
    ghost_cells(_ghost_cells)
{
    /* Each edge gains a collection of ghost cells for
     * boundary handling */
    total_cells = cells + 2 * ghost_cells;

    /* Set the constant step size */
    dx = calculate_dx();

    /* Initialise the field storage array */
    field = Vec1D(total_cells);

    /* Initialise and calculate the centroid coordinates storage array */
    cell_centres = Vec1D(total_cells);
    calculate_cell_centres();
}

/* Destructor */
CellCentredFieldMesh::~CellCentredFieldMesh()
{}

/* Obtain the field value at a particular cell index */
double CellCentredFieldMesh::operator() (const unsigned int index)
{
    return field[index];
}

/* Field value write access */
void CellCentredFieldMesh::set_field_value(const unsigned int index,
                                           const double _field_value)
{
    field[index] = _field_value;
}

/* Field value read access */
const double CellCentredFieldMesh::get_field_value(const unsigned int index) const
{
    return field[index];
}

/* Read access to the number of cells */
const unsigned int CellCentredFieldMesh::get_cells() const
{
    return cells;
}

/* Read access to the number of ghost cells on each edge */
const unsigned int CellCentredFieldMesh::get_ghost_cells() const
{
    return ghost_cells;
}

/* Read access to the total number of cells including ghost cells */
const unsigned int CellCentredFieldMesh::get_total_cells() const
{
    return total_cells;   
}

/* Read access to the x-direction step size */
const double CellCentredFieldMesh::get_dx() const
{
    return dx;
}

/* Read access to left side x-coordinate */
const double CellCentredFieldMesh::get_x_left() const
{
    return x_left;   
}

/* Read access to right side x-coordinate */
const double CellCentredFieldMesh::get_x_right() const
{
    return x_right;
}

/* Write access to entire field */
void CellCentredFieldMesh::set_field(const Vec1D& _field)
{
    field = _field;
}

/* Read access to entire field */
const Vec1D& CellCentredFieldMesh::get_field() const
{
    return field;
}

/* Read access to cell-centre coordinates storage array */
const Vec1D& CellCentredFieldMesh::get_cell_centres() const
{
    return cell_centres;
}

/* Read access to cell-centre coordinate at cell index */
const double CellCentredFieldMesh::get_cell_centre(const unsigned int index) const
{
    return cell_centres[index];
}

/* Backward difference */
double CellCentredFieldMesh::delta_back(const unsigned int index)
{
    return field[index] - field[index - 1];
}

/* Forward difference */
double CellCentredFieldMesh::delta_forward(const unsigned int index)
{
    return field[index + 1] - field[index];
}

/* Central difference */
double CellCentredFieldMesh::delta_central(const unsigned int index)
{
    return field[index + 1] - field[index - 1];
}

/* Initialise field values */
void CellCentredFieldMesh::initialise_field_values(const std::shared_ptr<InitialCondition>& init_con)
{
    for (unsigned int i=ghost_cells; i<total_cells - ghost_cells; i++) {
        double init_field_value = init_con->operator()(cell_centres[i]);
        set_field_value(i, init_field_value);
    }
}

/* Calculate the step size */
double CellCentredFieldMesh::calculate_dx()
{
    return (get_x_right() - get_x_left()) / cells;
}

/* Calculate the cell centre coordinates */
void CellCentredFieldMesh::calculate_cell_centres()
{
    double coord_start_x = get_x_left() - (static_cast<double>(ghost_cells) - 0.5) * dx;

    for (unsigned int i=0; i<total_cells; i++) {
        if ((i > ghost_cells - 1) && (i < total_cells - ghost_cells)) {
            cell_centres[i] = coord_start_x + static_cast<double>(i) * dx;
        } else {
            cell_centres[i] = 0.0;
        }
    }
}

}

#endif
