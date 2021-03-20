#ifndef __POLYGON_CPP
#define __POLYGON_CPP

#include "polygon.h"

namespace HyperFlow {

/* Constructor */
Polygon::Polygon()
{}

/* Constructor with vertices provided */
Polygon::Polygon(const Vec2D& _vertices) : vertices(_vertices)
{
    edges = calculate_edges();
    edge_areas = calculate_edge_areas();
    volume = calculate_volume();
    centroid = calculate_centroid();

    num_vertices = vertices.size();
    num_edges = edges.size();
}

/* Destructor */
Polygon::~Polygon()
{}

/* Get the number of vertices in the polygon */
unsigned int Polygon::get_num_vertices()
{
    return num_vertices;
}

/* Get the number of edges in the polygon */
unsigned int Polygon::get_num_edges()
{
    return num_edges;
}

/* Get the vertices in the polygon */
Vec2D& Polygon::get_vertices()
{
    return vertices;
}

/* Get the vector of edges of the polygon */
EdgeVec1D& Polygon::get_edges()
{
    return edges;
}

/* Get the edge 'areas' (lengths) */
Vec1D& Polygon::get_edge_areas()
{
    return edge_areas;
}

/* Get the polygon 'volume' (area) */
double Polygon::get_volume()
{
    return volume;
}

/* Get the centroid coordinates of the polygon */
Vec1D& Polygon::get_centroid()
{
    return centroid;
}

/* Obtain the minimum cartesian bounding box as a polygon */
Polygon Polygon::minimum_cartesian_bounding_box() const
{
    double min_left = 1e16;
    double max_right = -1e16;
    double min_bottom = 1e16;
    double max_top = -1e16;

    for (const auto& v : vertices) {
        if (v[0] <= min_left) {
            min_left = v[0];
        }
        if (v[0] >= max_right) {
            max_right = v[0];
        }
        if (v[1] <= min_bottom) {
            min_bottom = v[1];
        }
        if (v[1] >= max_top) {
            max_top = v[1];
        }
    }

    Vec1D bottom_left {min_left, min_bottom};
    Vec1D bottom_right {max_right, min_bottom};
    Vec1D top_right {max_right, max_top};
    Vec1D top_left {min_left, max_top};

    Vec2D bounding_box {bottom_left, bottom_right, top_right, top_left};
    Polygon bounding_poly = Polygon(bounding_box);
    return bounding_poly;
}

/* Output the polygon to the console */
void Polygon::output() const
{
    std::cout << std::setprecision(8) << "[\n";
    for (const auto& v : vertices) {
        v.output();
    }
    std::cout << "]" << std::endl;
}

/* Calculate the edges of the polygon */
EdgeVec1D Polygon::calculate_edges()
{
    EdgeVec1D edge_vec;

    for (unsigned int i=0; i<vertices.size(); i++) {
        if (i == vertices.size() - 1) {
            Edge e(vertices[i], vertices[0]);
            edge_vec.push_back(e);
        } else {
            Edge e(vertices[i], vertices[i+1]);
            edge_vec.push_back(e);
        }
    }

    return edge_vec;
}

/* Calculate the edge 'areas' (lengths) of the polygon */
Vec1D Polygon::calculate_edge_areas()
{
    Vec1D edge_lengths;

    for (const auto& edge : edges) {
        double edge_length = edge.get_length();
        edge_lengths.push_back(edge_length);
    }

    return edge_lengths;
}

/* Calculate the 'volume' of the polygon */
double Polygon::calculate_volume()
{
    double vol = 0.0;

    for (unsigned int i=0; i<vertices.size(); i++) {
        if (i == vertices.size() - 1) {
            vol += (vertices[i][0] * vertices[0][1] - vertices[0][0] * vertices[i][1]);
        } else {
            vol += (vertices[i][0] * vertices[i+1][1] - vertices[i+1][0] * vertices[i][1]);
        }
    }
    vol *= 0.5;

    return vol;
}

/* Calculate the centroid coordinates of the polygon */
Vec1D Polygon::calculate_centroid()
{
    double c_x = 0.0;
    double c_y = 0.0;

    for (unsigned int i=0; i<vertices.size(); i++) {
        double second = 0.0;
        if (i == vertices.size() - 1) {
            second = (vertices[i][0] * vertices[0][1] - vertices[0][0] * vertices[i][1]);
            c_x += ((vertices[i][0] + vertices[0][0]) * second);
            c_y += ((vertices[i][1] + vertices[0][1]) * second);
        } else {
            second = (vertices[i][0] * vertices[i+1][1] - vertices[i+1][0] * vertices[i][1]);
            c_x += ((vertices[i][0] + vertices[i+1][0]) * second);
            c_y += ((vertices[i][1] + vertices[i+1][1]) * second);
        }
    }

    c_x /= (6.0 * volume);
    c_y /= (6.0 * volume);

    Vec1D centroid_point {c_x, c_y};
    return centroid_point;
}

}

#endif
