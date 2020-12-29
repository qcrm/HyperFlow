#ifndef __EULER_DATA_OUTPUT_H
#define __EULER_DATA_OUTPUT_H

#include <fstream>

#include "../../../model/euler_ideal_gas/euler_ideal_gas.h"
#include "../data_output/data_output.h"

namespace HyperFlow {

class EulerIdealGasDataOutput
:
    public DataOutput
{
    
    public:
        
        /* Constructor */
        EulerIdealGasDataOutput();
        
        /* Constructor supplying spatio-temporal extents
         * for results metadata */
        EulerIdealGasDataOutput(const std::shared_ptr<Model> _model,
                                const double _x_left,
                                const double _x_right,
                                const double _y_bottom,
                                const double _y_top,
                                const double _t_start,
                                const double _t_end);
        
        /* Destructor */
        virtual ~EulerIdealGasDataOutput();
        
        /* Output data to file */
        virtual void to_file(const std::shared_ptr<Mesh>& mesh,
                     const unsigned int n,
                     const double t);
};

}

#endif
