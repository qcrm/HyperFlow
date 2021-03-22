#ifndef __GODUNOV_CPP
#define __GODUNOV_CPP

#include "godunov.h"
#include <iostream>

namespace HyperFlow {

/* Constructor */
GodunovScheme::GodunovScheme()
{}

/* Parameterised constructor */
GodunovScheme::GodunovScheme(
    const std::shared_ptr<Model> _model,
    const std::shared_ptr<RiemannSolver> _riemann,
    std::shared_ptr<Mesh> _mesh)
:
    model(_model),
    riemann(_riemann),
    mesh(_mesh)
{}

/* Destructor */
GodunovScheme::~GodunovScheme()
{}

/* Carry out the Godunov scheme on the provided mesh */
Vec2D GodunovScheme::operator()()
{
    Vec2D cell_updates = Vec2D(mesh->get_cells().size(), Vec1D(4, 0.0));

    /* Rotated flux calculations */
    for (unsigned int cell_edge_idx=0; cell_edge_idx<mesh->get_cell_edge_indices().size(); cell_edge_idx++) {
        int cell_idx = mesh->get_cell_edge_indices()[cell_edge_idx][0];
        int nbrc_idx = mesh->get_cell_edge_indices()[cell_edge_idx][1];
        int edge_idx = mesh->get_cell_edge_indices()[cell_edge_idx][2];

        Cell cell = mesh->get_cells()[cell_idx];
        
        Cell nbrc;
        if (nbrc_idx == -1) {  // Solid
            nbrc = cell;
        } else if (nbrc_idx == -2) {  // Transmissive
            nbrc = cell;
        } else {
            nbrc = mesh->get_cells()[nbrc_idx];
        }
        Edge edge = mesh->get_edges()[edge_idx];

        double normal_angle = edge.get_outward_normal_angle();

        Vec1D cell_values = cell.get_flow_values();
        Vec1D nbrc_values = nbrc.get_flow_values();

        Vec1D rot_cell_values = model->rotate_flow_values(normal_angle, cell_values);
        Vec1D rot_nbrc_values = model->rotate_flow_values(normal_angle, nbrc_values);

        if (nbrc_idx == -1) {
            rot_nbrc_values[1] *= -1.0;
        }

        Vec1D rot_flux = riemann->operator()(rot_cell_values, rot_nbrc_values, Direction::x);
        Vec1D flux = model->rotate_flow_values(-1.0 * normal_angle, rot_flux);
        flux *= (edge.get_length());
        cell_updates[cell_idx] += flux;
    }

    /* Cell spatial updates */
    for (unsigned int cell_idx=0; cell_idx<mesh->get_cells().size(); cell_idx++) {
        Cell cell = mesh->get_cells()[cell_idx];

        double inv_cell_volume = 1.0 / cell.get_volume();
        
        cell_updates[cell_idx] *= inv_cell_volume;
        cell_updates[cell_idx] *= -1.0;
    }

    return cell_updates;
}

}

#endif
