#ifndef __CELL_H
#define __CELL_H

#include "../../geometry/polygon/polygon.h"

namespace HyperFlow {

class Cell
:
    public Polygon
{
    
public:

    /* Constructor */
    Cell();
    
    /* Constructor with vertices */
    Cell(const Vec2D& _vertices);
    
    /* Constructor with vertices and flow values */
    Cell(const Vec2D& _vertices, const Vec1D& _flow_values);
    
    /* Destructor */
    virtual ~Cell();
    
    /* Get the flow values vector */
    const Vec1D& get_flow_values() const;
    
    /* Set the flow values vector */
    void set_flow_values(const Vec1D& flow_values);

protected:
    
    /* The vector of flow values */
    Vec1D flow_values;

};

typedef std::vector<Cell> CellVec1D;

}

#endif
