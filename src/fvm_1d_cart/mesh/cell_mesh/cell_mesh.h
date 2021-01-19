#ifndef __CELL_CENTRED_FIELD_MESH_H
#define __CELL_CENTRED_FIELD_MESH_H

#include "../../tensor/tensor.h"
#include "../../simulation/initcon/initcon/initcon.h"

namespace HyperFlow {

class CellCentredFieldMesh {

    public:

        /* Constructor */
        CellCentredFieldMesh();

        /* Constructor from number of cells */
        CellCentredFieldMesh(const double _x_left,
                             const double _x_right,
                             const unsigned int _cells,
                             const unsigned int _ghost_cells);

        /* Destructor */
        virtual ~CellCentredFieldMesh();

        /* Obtain the field value at a particular cell index */
        double operator() (const unsigned int index);

        /* Field value write access */
        void set_field_value(const unsigned int index,
                             const double _field_value);

        /* Field value read access */
        const double get_field_value(const unsigned int index) const;

        /* Read access to the number of cells */
        const unsigned int get_cells() const;

        /* Read access to the number of ghost cells on each edge */
        const unsigned int get_ghost_cells() const;

        /* Read access to the total number of cells including ghost cells */
        const unsigned int get_total_cells() const;

        /* Read access to the x-direction step size */
        const double get_dx() const;

        /* Read access to left side x-coordinate */
        const double get_x_left() const;

        /* Read access to right side x-coordinate */
        const double get_x_right() const;

        /* Write access to entire field */
        void set_field(const Vec1D& _field);

        /* Read access to entire field */
        const Vec1D& get_field() const;

        /* Read access to cell-centre coordinates storage array */
        const Vec1D& get_cell_centres() const;

        /* Read access to cell-centre coordinate at cell index */
        const double get_cell_centre(const unsigned int index) const;

        /* Backward difference */
        double delta_back(const unsigned int index);

        /* Forward difference */
        double delta_forward(const unsigned int index);

        /* Central difference */
        double delta_central(const unsigned int index);

        /* Initialise field values */
        void initialise_field_values(const std::shared_ptr<InitialCondition>& init_con);

    private:

        /* Left x coordinate */
        double x_left;

        /* Right x coordinate */
        double x_right;

        /* Number of non ghost cells */
        unsigned int cells;

        /* Number of ghost cells for each edge */
        unsigned int ghost_cells;

        /* Total cells (including ghost cells) */
        unsigned int total_cells;

        /* Homogeneous step-size */
        double dx;

        /* Array storing field values */
        Vec1D field;

        /* Coordinates of cell centres */
        Vec1D cell_centres;

        /* Calculate the step size */
        double calculate_dx();

        /* Calculate the cell centre coordinates */
        void calculate_cell_centres();

};

}

#endif
