#ifndef __OUTPUT_CPP
#define __OUTPUT_CPP

#include "output.h"

namespace HyperFlow {

/* Constructor */
EulerIdealGasDataOutput::EulerIdealGasDataOutput()
{}

/* Constructor supplying spatio-temporal extents
 * for results metadata */
EulerIdealGasDataOutput::EulerIdealGasDataOutput(const std::shared_ptr<EulerIdealGasModel> _model,
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
EulerIdealGasDataOutput::~EulerIdealGasDataOutput() {};

/* Output data to file */
void EulerIdealGasDataOutput::to_file(const std::shared_ptr<Mesh>& mesh,
                                      const unsigned int n,
                                      const double t)
{
    char filename[19];
    sprintf(filename, "results_%06d.csv", n);
    std::string header = "x,y,rho,rho_u,rho_v,E,u,v,p,e,M";
 
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

    std::vector<std::shared_ptr<MeshBlock> > mesh_blocks = mesh->get_mesh_blocks();

    unsigned int ghost_cells = mesh_blocks[0]->get_ghost_cells();

    for (unsigned int b=0; b<mesh_blocks.size(); b++) {
        for (unsigned int i=ghost_cells; i<mesh_blocks[b]->get_total_x_cells() - ghost_cells; i++) {
            for (unsigned int j=ghost_cells; j<mesh_blocks[b]->get_total_y_cells() - ghost_cells; j++) {                
                Vec1D centroid = mesh_blocks[b]->get_centroid(i, j);
                Vec1D cons = mesh_blocks[b]->get_flow_values(i, j);
                Vec1D prim = model->cons_to_prim(cons);

                double mach_number = model->prim_mach_number(prim);

                out << centroid[0] << ",";         // x-coordinate
                out << centroid[1] << ",";         // y-coordinate
                out << cons[0] << ",";                            // Density (rho)
                out << cons[1] << ",";                            // x-Momentum (rho_u)
                out << cons[2] << ",";                            // y-Momentum (rho_v)
                out << cons[3] << ",";                            // Total energy (E)
                out << prim[1] << ",";                            // x-Velocity (u)
                out << prim[2] << ",";                            // v-Velocity (v)
                out << prim[3] << ",";                            // Pressure (p)
                out << model->prim_internal_energy(prim) << ",";       // Internal energy (e)
                out << mach_number;                        // Mach number
                out << std::endl;
            }
        }
    }

    out.close();
}

}

#endif
