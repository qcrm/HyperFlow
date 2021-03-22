#ifndef __OUTPUT_H
#define __OUTPUT_H

#include <fstream>

#include "../../model/euler_ideal_gas/euler_ideal_gas.h"
#include "../../mesh/cell/cell.h"
#include "../../mesh/mesh/mesh.h"

namespace HyperFlow {

class EulerIdealGasDataOutput
{
    
    public:
        
        /* Constructor */
        EulerIdealGasDataOutput();
        
        /* Constructor supplying spatio-temporal extents
         * for results metadata */
        EulerIdealGasDataOutput(const std::shared_ptr<EulerIdealGasModel> _model,
                                const double _x_left,
                                const double _x_right,
                                const double _y_bottom,
                                const double _y_top,
                                const double _t_start,
                                const double _t_end);
        
        /* Destructor */
        virtual ~EulerIdealGasDataOutput();
        
        /* Output data to file */
        void to_file(const std::shared_ptr<Mesh>& mesh,
                     const unsigned int n,
                     const double t);

    private:
       
        /* Pointer to the model equations */
        std::shared_ptr<EulerIdealGasModel> model;
        
        /* Left domain extent */
        double x_left;
        
        /* Right domain extent */
        double x_right;
        
        /* Bottom domain extent */
        double y_bottom;
        
        /* Top domain extent */
        double y_top;
        
        /* Start time of the simulation */
        double t_start;
        
        /* End time of the simulation */
        double t_end;
};

}

#endif
