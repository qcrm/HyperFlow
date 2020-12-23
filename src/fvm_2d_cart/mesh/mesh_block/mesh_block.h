#ifndef __MESH_BLOCK_H
#define __MESH_BLOCK_H

#include "../../tensor/tensor.h"
#include "../../simulation/initcon/initcon/initcon.h"

namespace HyperFlow {

class MeshBlock {

    public:

        /* Constructor */
        MeshBlock();

        /* Constructor from number of cells */
        MeshBlock(const Vec1D& _extents,
                  const unsigned int _x_cells,
                  const unsigned int _y_cells,
                  const unsigned int _ghost_cells,
                  const unsigned int _dimension);

        /* Destructor */
        virtual ~MeshBlock();

        /* Obtain a vector of the field values at a particular index */
        Vec1D& operator() (const unsigned int x_index,
                           const unsigned int y_index);

        /* Flow value write access */
        void set_flow_values(const unsigned int x_index,
                             const unsigned int y_index,
                             const Vec1D& _flow_values);

        /* Flow value read access */
        const Vec1D& get_flow_values(const unsigned int x_index,
                                     const unsigned int y_index) const;

        /* Read access to the number of x cells */
        const unsigned int get_x_cells() const;

        /* Read access to the number of y cells */
        const unsigned int get_y_cells() const;

        /* Read access to the number of ghost cells on each edge */
        const unsigned int get_ghost_cells() const;

        /* Read access to the total number of cells in the x direction */
        const unsigned int get_total_x_cells() const;

        /* Read access to the total number of cells in the y direction */
        const unsigned int get_total_y_cells() const;

        /* Read access to the dimension of the model system */
        const unsigned int get_dimension() const;

        /* Read access to the x-direction step size */
        const double get_dx() const;

        /* Read access to the y-direction step size */
        const double get_dy() const;

        /* Read access to left side x-coordinate */
        const double get_x_left() const;

        /* Read access to right side x-coordinate */
        const double get_x_right() const;

        /* Read access to bottom side y-coordinate */
        const double get_y_bottom() const;

        /* Read access to top side y-coordinate */
        const double get_y_top() const;

        /* Write access to entire flow fields */
        void set_fields(const Vec3D& _fields);

        /* Read access to entire flow fields */
        const Vec3D& get_fields() const;

        /* Read access to centroid coordinates storage array */
        const Vec3D& get_centroids() const;

        /* Read access to centroid coordinates at cell index */
        const Vec1D& get_centroid(const unsigned int x_index,
                                  const unsigned int y_index) const;

        /* Backward difference in the x-direction */
        Vec1D delta_x_back(const unsigned int x_index,
                           const unsigned int y_index);

        /* Forward difference in the x-direction */
        Vec1D delta_x_forward(const unsigned int x_index,
                              const unsigned int y_index);


        /* Central difference in the x-direction */
        Vec1D delta_x_central(const unsigned int x_index,
                              const unsigned int y_index);


        /* Backward difference in the y-direction */
        Vec1D delta_y_back(const unsigned int x_index,
                           const unsigned int y_index);


        /* Forward difference in the y-direction */
        Vec1D delta_y_forward(const unsigned int x_index,
                              const unsigned int y_index);


        /* Central difference in the y-direction */
        Vec1D delta_y_central(const unsigned int x_index,
                              const unsigned int y_index);

        /* Initialise field values */
        void initialise_field_values(const std::shared_ptr<InitialCondition>& init_con);

    private:

        /* Vector containing extents of the rectangular domain:
         * [0] - x_left, [1] - x_right, [2] - y_bottom, [3] - y_top */
        Vec1D extents;

        /* Number of non ghost cells in x-direction */
        unsigned int x_cells;

        /* Number of non ghost cells in y-direction */
        unsigned int y_cells;

        /* Number of ghost cells for each edge */
        unsigned int ghost_cells;

        /* Dimension of the model system */
        unsigned int dimension;

        /* Total cells (including ghost cells) in the x-direction */
        unsigned int total_x_cells;

        /* Total cells (including ghost cells) in the y-direciton */
        unsigned int total_y_cells;

        /* Homogeneous step-size in x-direction */
        double dx;

        /* Homogeneous step-size in y-direction */
        double dy;

        /* Multidimensional array storing flow values for mesh block */
        Vec3D fields;

        /* Coordinates of cell centroids */
        Vec3D centroids;

        /* Calculate the x-direction step size */
        double calculate_dx();

        /* Calculate the y-direction step size */
        double calculate_dy();

        /* Calculate the cell centroid coordinates */
        void calculate_centroids();

};

}

#endif
