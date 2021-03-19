#ifndef __EDGE_CPP
#define __EDGE_CPP

#include "edge.h"

namespace HyperFlow {

/* Constructor */
Edge::Edge()
{}

/* Constructor with start, end points and distance tolerance */
Edge::Edge(const Vec1D& _start,
           const Vec1D& _end,
           double _tol)
:
    start(_start),
    end(_end),
    tol(_tol)
{
    midpoint = calculate_midpoint();
    point_vector = calculate_point_vector();
    inward_unit_normal = calculate_inward_unit_normal();
    outward_unit_normal = calculate_outward_unit_normal();
    point_vector_angle = calculate_point_vector_angle();
    inward_normal_angle = calculate_inward_normal_angle();
    outward_normal_angle = calculate_outward_normal_angle();
    horizontal = calculate_horizontal();
    vertical = calculate_vertical();
    axis_aligned = calculate_axis_aligned();
}

/* Destructor */
Edge::~Edge()
{}

/* Get starting point of edge */ 
const Vec1D& Edge::get_start() const
{
    return start;
}

/* Get ending point of edge */ 
const Vec1D& Edge::get_end() const
{
    return end;
}

/* Get the midpoint coordinates of the edge */
const Vec1D& Edge::get_midpoint() const
{
    return midpoint;
}

/* Get the length of the edge */
const double Edge::get_length() const
{
    return point_vector.l2_norm();
}

/* Get the vector from the start to the end of the edge */
const Vec1D& Edge::get_point_vector()
{
    return point_vector;
}

/* Get the inward unit normal vector */
const Vec1D& Edge::get_inward_unit_normal()
{
    return inward_unit_normal;
}

/* Get the outward unit normal vector */
const Vec1D& Edge::get_outward_unit_normal()
{
    return outward_unit_normal;
}

/* Get the angle of the 'point' vector */
const double Edge::get_point_vector_angle()
{
    return point_vector_angle;
}

/* Get the inward normal angle */
const double Edge::get_inward_normal_angle()
{
    return inward_normal_angle;
}

/* Get the outward normal angle */
const double Edge::get_outward_normal_angle()
{
    return outward_normal_angle;
}

/* Is the edge vertical? That is, aligned with the y coordinate axis? */
const bool Edge::is_vertical()
{
    return vertical;
}

/* Is the edge horizontal? That is, aligned with the x coordinate axis? */
const bool Edge::is_horizontal()
{
    return horizontal;
}

/* Is the edge aligned to either of the coordinate axes? */
const bool Edge::is_axis_aligned()
{
    return axis_aligned;
}

/* Calculate the midpoint coordinates of the edge */ 
Vec1D Edge::calculate_midpoint()
{
    Vec1D midpoint {(end.x + start.x) / 2.0, (end.y + start.y) / 2.0};
    return midpoint;
}

/* Calculate the End - Start vector of the edge */
Vec1D Edge::calculate_point_vector()
{
    Vec1D p = end - start;
    Vec1D v(2);
    v[0] = p.x;
    v[1] = p.y;
    return v;
}

/* Calculate the inward unit normal of the edge */
Vec1D Edge::calculate_inward_unit_normal()
{
    Vec1D normal(2);
    normal[0] = -1.0 * point_vector[1];
    normal[1] = point_vector[0];

    double normal_length = sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
    normal[0] /= normal_length;
    normal[1] /= normal_length;

    return normal;
}

/* Calculate the outward unit normal of the edge */
Vec1D Edge::calculate_outward_unit_normal()
{
    Vec1D normal(2);
    normal[0] = point_vector[1];
    normal[1] = -1.0 * point_vector[0];

    double normal_length = sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
    normal[0] /= normal_length;
    normal[1] /= normal_length;

    return normal;
}

/* Calculate the End - Start vector angle */
double Edge::calculate_point_vector_angle()
{
    return atan2(point_vector[1], point_vector[0]) * 180.0 / M_PI;
}

/* Calculate the inward unit normal vector angle */ 
double Edge::calculate_inward_normal_angle()
{
    return atan2(inward_unit_normal[1], inward_unit_normal[0]) * 180.0 / M_PI;
}

/* Calculate the outward unit normal vector angle */
double Edge::calculate_outward_normal_angle()
{
    return atan2(outward_unit_normal[1], outward_unit_normal[0]) * 180.0 / M_PI;
}

/* Calculate whether the edge is horizontal */
bool Edge::calculate_horizontal()
{
    return ((fabs(point_vector[0]) > tol) && (fabs(point_vector[1]) < tol));
}

/* Calculate whether the edge is vertical */
bool Edge::calculate_vertical()
{
    return ((fabs(point_vector[0]) < tol) && (fabs(point_vector[1]) > tol));
}

/* Calculate whether the edge is axis-aligned */
bool Edge::calculate_axis_aligned()
{
    return (vertical || horizontal);
}

}

#endif
