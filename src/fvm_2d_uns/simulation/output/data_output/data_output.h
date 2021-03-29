#ifndef __DATA_OUTPUT_H
#define __DATA_OUTPUT_H

#include <fstream>

#include "../../../model/model/model.h"
#include "../../../mesh/mesh/mesh.h"

namespace HyperFlow {

class DataOutput
{
    
    public:
        
        /* Constructor */
        DataOutput();
        
        /* Constructor supplying spatio-temporal extents
         * for results metadata */
        DataOutput(const std::shared_ptr<Model> _model,
                   const double _x_left,
                   const double _x_right,
                   const double _y_bottom,
                   const double _y_top,
                   const double _t_start,
                   const double _t_end);
        
        /* Destructor */
        virtual ~DataOutput();
        
        /* Output data to file */
        virtual void to_file(const std::shared_ptr<Mesh>& mesh,
                             const unsigned int n,
                             const double t) = 0;

    protected:
       
        /* Pointer to the model equations */
        std::shared_ptr<Model> model;
        
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
