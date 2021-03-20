#ifndef __POLYGON_H
#define __POLYGON_H

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>

#include "../edge/edge.h"

namespace HyperFlow {

class Polygoni
{

public:

    /* Constructor */
    Polygon();

    /* Constructor with vertices provided */
    Polygon(const Vec2D& _vertices);

    /* Destructor */
    virtual ~Polygon();

    /* Get the number of vertices in the polygon */
    unsigned int get_num_vertices();
    
    /* Get the number of edges in the polygon */
    unsigned int get_num_edges();

    /* Get the vertices in the polygon */
    Vec2D& get_vertices();

    /* Get the vector of edges of the polygon */
    EdgeVec1D& get_edges();

    /* Get the edge 'areas' (lengths) */
    Vec1D& get_edge_areas();
    
    /* Get the polygon 'volume' (area) */
    double get_volume();

    /* Get the centroid coordinates of the polygon */
    Vec1D& get_centroid();

    /* Obtain the minimum cartesian bounding box as a polygon */
    Polygon minimum_cartesian_bounding_box() const;

    /* Output the polygon to the console */
    void output() const;

protected:

    /* Number of vertices in the polygon */
    unsigned int num_vertices;

    /* Number of edges in the polygon */
    unsigned int num_edges;

    /* The vertex coordinates list for the polygon */
    Vec2D vertices;
    
    /* The vector of edges for the polygon */
    EdgeVec1D edges;

    /* The edge 'areas' (lengths in 2D) */
    Vec1D edge_areas;

    /* The polygon 'volume' (area in 2D) */
    double volume;

    /* The centroid coordinates of the polygon */
    Vec1D centroid;
    
    /* Calculate the edges of the polygon */
    EdgeVec1D calculate_edges();

    /* Calculate the edge 'areas' (lengths) of the polygon */
    Vec1D calculate_edge_areas();

    /* Calculate the 'volume' of the polygon */
    double calculate_volume();

    /* Calculate the centroid coordinates of the polygon */
    Vec1D calculate_centroid();

};

#endif
