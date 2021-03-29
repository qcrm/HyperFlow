#ifndef __DATA_OUTPUT_CPP
#define __DATA_OUTPUT_CPP

#include "data_output.h"

namespace HyperFlow {

/* Constructor */
DataOutput::DataOutput()
{}

/* Constructor supplying spatio-temporal extents
 * for results metadata */
DataOutput::DataOutput(const std::shared_ptr<Model> _model,
                       const double _x_left,
                       const double _x_right,
                       const double _y_bottom,
                       const double _y_top,
                       const double _t_start,
                       const double _t_end)
:
    model(_model),
    x_left(_x_left),
    x_right(_x_right),
    y_bottom(_y_bottom),
    y_top(_y_top),
    t_start(_t_start),
    t_end(_t_end)
{}

/* Destructor */
DataOutput::~DataOutput() {};

}

#endif
