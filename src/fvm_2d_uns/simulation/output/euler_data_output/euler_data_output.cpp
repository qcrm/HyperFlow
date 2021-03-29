#ifndef __EULER_DATA_OUTPUT_CPP
#define __EULER_DATA_OUTPUT_CPP

#include "euler_data_output.h"

namespace HyperFlow {

/* Constructor */
EulerIdealGasDataOutput::EulerIdealGasDataOutput()
{}

/* Constructor supplying spatio-temporal extents
 * for results metadata */
EulerIdealGasDataOutput::EulerIdealGasDataOutput(const std::shared_ptr<Model> _model,
                                                 const double _x_left,
                                                 const double _x_right,
                                                 const double _y_bottom,
                                                 const double _y_top,
                                                 const double _t_start,
                                                 const double _t_end)
:
    DataOutput(
        _model,
        _x_left,
        _x_right,
        _y_bottom,
        _y_top,
        _t_start,
        _t_end
    )
{}

/* Destructor */
EulerIdealGasDataOutput::~EulerIdealGasDataOutput() {};

/* Output data to file */
void EulerIdealGasDataOutput::to_file(const std::shared_ptr<Mesh>& mesh,
                                      const unsigned int n,
                                      const double t)
{
    char filename[19];
    sprintf(filename, "results_%06d.csv", n);
    std::string header = "x,y,rho,rho_u,rho_v,E,u,v,p";
 
    std::ofstream out(filename);

    // Add metadata
    out << "X_LEFT=" << x_left << std::endl;
    out << "X_RIGHT=" << x_right << std::endl;
    out << "Y_BOTTOM=" << y_bottom << std::endl;
    out << "Y_TOP=" << y_top << std::endl;
    out << "T_START=" << t_start << std::endl;
    out << "T_END=" << t_end << std::endl;
    out << "T=" << t << std::endl << std::endl;

    out << header << std::endl; 

    CellVec1D cells = mesh->get_cells();

    for (unsigned int cell_idx=0; cell_idx<cells.size(); cell_idx++) {
        Cell cell = cells[cell_idx];
        Vec1D centroid = cell.get_centroid();
        Vec1D cons = cell.get_flow_values();
        Vec1D prim = model->cons_to_prim(cons);

        out << centroid[0] << ",";                          // x-coordinate
        out << centroid[1] << ",";                          // y-coordinate
        out << cons[0] << ",";                              // Density (rho)
        out << cons[1] << ",";                              // x-Momentum (rho_u)
        out << cons[2] << ",";                              // y-Momentum (rho_v)
        out << cons[3] << ",";                              // Total energy (E)
        out << prim[1] << ",";                              // x-Velocity (u)
        out << prim[2] << ",";                              // v-Velocity (v)
        out << prim[3];                                     // Pressure (p)
        out << std::endl;
    }

    out.close();
}

}

#endif
