#ifndef __MESH_BLOCK_CPP
#define __MESH_BLOCK_CPP

#include "mesh_block.h"

namespace HyperFlow {

/* Constructor */
MeshBlock::MeshBlock()
{}

/* Constructor from number of cells */
MeshBlock::MeshBlock(const Vec1D& _extents,
                     const unsigned int _x_cells,
                     const unsigned int _y_cells,
                     const unsigned int _ghost_cells,
                     const unsigned int _dimension)
:
    extents(_extents),
    x_cells(_x_cells),
    y_cells(_y_cells),
    ghost_cells(_ghost_cells),
    dimension(_dimension)
{
    /* Each edge gains a collection of ghost cells for
     * boundary handling */
    total_x_cells = x_cells + 2 * ghost_cells;
    total_y_cells = y_cells + 2 * ghost_cells;

    /* Set the constant step sizes in each direction */
    dx = calculate_dx();
    dy = calculate_dy();

    /* Initialise the fields storage array */
    fields = Vec3D(total_x_cells, Vec2D(total_y_cells, Vec1D(dimension)));

    /* Initialise and calculate the centroid coordinates storage array */
    centroids = Vec3D(total_x_cells, Vec2D(total_y_cells, Vec1D(2)));
    calculate_centroids();
}

/* Destructor */
MeshBlock::~MeshBlock()
{}

/* Obtain a vector of the field values at a particular index */
Vec1D& MeshBlock::operator() (const unsigned int x_index,
                              const unsigned int y_index)
{
    return fields[x_index][y_index];
}


/* Flow value write access */
void MeshBlock::set_flow_values(const unsigned int x_index,
                                const unsigned int y_index,
                                const Vec1D& _flow_values)
{
    fields[x_index][y_index] = _flow_values;
}

/* Flow value read access */
const Vec1D& MeshBlock::get_flow_values(const unsigned int x_index,
                                        const unsigned int y_index) const
{
    return fields[x_index][y_index];
}

/* Read access to the number of x cells */
const unsigned int MeshBlock::get_x_cells() const
{
    return x_cells;
}

/* Read access to the number of y cells */
const unsigned int MeshBlock::get_y_cells() const
{
    return y_cells;
}

/* Read access to the number of ghost cells on each edge */
const unsigned int MeshBlock::get_ghost_cells() const
{
    return ghost_cells;
}

/* Read access to the total number of cells in the x direction */
const unsigned int MeshBlock::get_total_x_cells() const
{
    return total_x_cells;
}

/* Read access to the total number of cells in the y direction */
const unsigned int MeshBlock::get_total_y_cells() const
{
    return total_y_cells;
}

/* Read access to the dimension of the model system */
const unsigned int MeshBlock::get_dimension() const
{
    return dimension;
}

/* Read access to the x-direction step size */
const double MeshBlock::get_dx() const
{
    return dx;
}

/* Read access to the y-direction step size */
const double MeshBlock::get_dy() const
{
    return dy;
}

/* Read access to left side x-coordinate */
const double MeshBlock::get_x_left() const
{
    return extents[0];
}

/* Read access to right side x-coordinate */
const double MeshBlock::get_x_right() const
{
    return extents[1];
}

/* Read access to bottom side y-coordinate */
const double MeshBlock::get_y_bottom() const
{
    return extents[2];
}

/* Read access to top side y-coordinate */
const double MeshBlock::get_y_top() const
{
    return extents[3];
}

/* Read access to centroids */
const Vec3D& MeshBlock::get_centroids() const
{
    return centroids;
}

/* Calculation of centroid for a cell index */
const Vec1D& MeshBlock::get_centroid(const unsigned int x_index,
                                     const unsigned int y_index) const
{
    return centroids[x_index][y_index];
}

/* Write access to entire flow fields */
void MeshBlock::set_fields(const Vec3D& _fields)
{
    fields = _fields;
}

/* Read access to entire flow fields */
const Vec3D& MeshBlock::get_fields() const
{
    return fields;
}

/* Backward difference in the x-direction */
Vec1D MeshBlock::delta_x_back(const unsigned int x_index,
                              const unsigned int y_index)
{
    Vec1D delta = fields[x_index][y_index] - fields[x_index - 1][y_index];
    return delta;
}

/* Forward difference in the x-direction */
Vec1D MeshBlock::delta_x_forward(const unsigned int x_index,
                                 const unsigned int y_index)
{
    Vec1D delta = fields[x_index + 1][y_index] - fields[x_index][y_index];
    return delta;
}

/* Central difference in the x-direction */
Vec1D MeshBlock::delta_x_central(const unsigned int x_index,
                                 const unsigned int y_index)
{
    Vec1D delta = fields[x_index + 1][y_index] - fields[x_index - 1][y_index];
    return delta;
}

/* Backward difference in the y-direction */
Vec1D MeshBlock::delta_y_back(const unsigned int x_index,
                              const unsigned int y_index)
{
    Vec1D delta = fields[x_index][y_index] - fields[x_index][y_index - 1];
    return delta;
}

/* Forward difference in the y-direction */
Vec1D MeshBlock::delta_y_forward(const unsigned int x_index,
                                 const unsigned int y_index)
{
    Vec1D delta = fields[x_index][y_index + 1] - fields[x_index][y_index];
    return delta;
}

/* Central difference in the y-direction */
Vec1D MeshBlock::delta_y_central(const unsigned int x_index,
                                 const unsigned int y_index)
{
    Vec1D delta = fields[x_index][y_index + 1] - fields[x_index][y_index - 1];
    return delta;
}

/* Initialise field values */
void MeshBlock::initialise_field_values(const std::shared_ptr<InitialCondition>& init_con)
{
    for (unsigned int i=ghost_cells; i<total_x_cells - ghost_cells; i++) {
        for (unsigned int j=ghost_cells; j<total_y_cells - ghost_cells; j++) {
            Vec1D init_flow_values = init_con->operator()(centroids[i][j][0],
                                                          centroids[i][j][1]);
            set_flow_values(i, j, init_flow_values);;
        }
    }
}

/* Calculate the x-direction step size */
double MeshBlock::calculate_dx()
{
    return (get_x_right() - get_x_left()) / x_cells;
}

/* Calculate the y-direction step size */
double MeshBlock::calculate_dy()
{
    return (get_y_top() - get_y_bottom()) / y_cells;
}

/* Calculate the cell centroid coordinates */
void MeshBlock::calculate_centroids() {
    double coord_start_x = get_x_left() - (static_cast<double>(ghost_cells) - 0.5) * dx;
    double coord_start_y = get_y_bottom() - (static_cast<double>(ghost_cells) - 0.5) * dy;

    for (unsigned int i=0; i<total_x_cells; i++) {
        for (unsigned int j=0; j<total_y_cells; j++) {
            if (
                (i > ghost_cells - 1) && (i < total_x_cells - ghost_cells) &&
                (j > ghost_cells - 1) && (j < total_y_cells - ghost_cells)
            ) {
                centroids[i][j][0] = coord_start_x + static_cast<double>(i) * dx;
                centroids[i][j][1] = coord_start_y + static_cast<double>(j) * dy;
            } else {
                centroids[i][j][0] = 0.0;
                centroids[i][j][1] = 0.0;
            }
        }
    }
}

}

#endif
