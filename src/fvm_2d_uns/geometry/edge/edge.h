#ifndef __EDGE_H
#define __EDGE_H

#include <cmath>

#include "../../../share/tensor/tensor.h"

namespace HyperFlow {

class Edge
{

public:
    
    /* Constructor */
    Edge();

    /* Constructor with start, end points and distance tolerance */
    Edge(const Vec1D& _start,
         const Vec1D& _end,
         const double _tol = 1e-12);

    /* Destructor */
    virtual ~Edge();

    /* Are two edges equal? */
    bool operator==(const Edge& rhs);

    /* Are two edges equal in 'opposite' direction? */
    bool opposite_edges_equal(const Edge& rhs);

    /* Get starting point of edge */
    const Vec1D& get_start() const;

    /* Get ending point of edge */
    const Vec1D& get_end() const;

    /* Get the midpoint coordinates of the edge */
    const Vec1D& get_midpoint() const;

    /* Get the length of the edge */
    const double get_length() const;

    /* Get the vector from the start to the end of the edge */
    const Vec1D& get_point_vector();

    /* Get the inward unit normal vector */
    const Vec1D& get_inward_unit_normal();

    /* Get the outward unit normal vector */
    const Vec1D& get_outward_unit_normal();

    /* Get the angle of the 'point' vector */
    const double get_point_vector_angle();

    /* Get the inward normal angle */
    const double get_inward_normal_angle();

    /* Get the outward normal angle */
    const double get_outward_normal_angle();

    /* Is the edge vertical? That is, aligned with the y coordinate axis? */
    const bool is_vertical();

    /* Is the edge horizontal? That is, aligned with the x coordinate axis? */
    const bool is_horizontal();

    /* Is the edge aligned to either of the coordinate axes? */
    const bool is_axis_aligned();

private:

    /* Starting coordinates of the edge */
    Vec1D start;

    /* Ending coordinates of the edge */
    Vec1D end;

    /* Midpoint coordinates of the edge */
    Vec1D midpoint;

    /* Tolerance value used to control equality */
    double tol;

    /* End - Start vector */
    Vec1D point_vector;

    /* Inward unit normal vector of the edge */
    Vec1D inward_unit_normal;

    /* Outward unit normal vector of the edge */
    Vec1D outward_unit_normal;

    /* Angle of the edge w.r.t. horizontal coordinate x axis */ 
    double point_vector_angle;

    /* Angle of the inward normal w.r.t. horizontal coordinate x axis */ 
    double inward_normal_angle;

    /* Angle of the outward normal w.r.t. horizontal coordinate x axis */ 
    double outward_normal_angle;

    /* Is the edge vertical? */
    bool vertical;

    /* Is the edge horizontal? */
    bool horizontal;

    /* Is the edge axis-aligned? */
    bool axis_aligned;

    /* Calculate the midpoint coordinates of the edge */
    Vec1D calculate_midpoint();

    /* Calculate the End - Start vector of the edge */
    Vec1D calculate_point_vector();

    /* Calculate the inward unit normal of the edge */
    Vec1D calculate_inward_unit_normal();

    /* Calculate the outward unit normal of the edge */
    Vec1D calculate_outward_unit_normal();

    /* Calculate the End - Start vector angle */
    double calculate_point_vector_angle();

    /* Calculate the inward unit normal vector angle */
    double calculate_inward_normal_angle();

    /* Calculate the outward unit normal vector angle */
    double calculate_outward_normal_angle();

    /* Calculate whether the edge is horizontal */
    bool calculate_horizontal();

    /* Calculate whether the edge is vertical */
    bool calculate_vertical();

    /* Calculate whether the edge is axis-aligned */
    bool calculate_axis_aligned();

};

typedef std::vector<Edge> EdgeVec1D;

}

#endif
