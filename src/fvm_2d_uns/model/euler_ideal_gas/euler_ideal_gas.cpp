#ifndef __EULER_IDEAL_GAS_CPP
#define __EULER_IDEAL_GAS_CPP

#include "euler_ideal_gas.h"

namespace HyperFlow {

/* Constructor */
EulerIdealGasModel::EulerIdealGasModel()
{}

/* Constructor via ratio of specific heats,
 * defaulting to air (\gamma = 1.4) */
EulerIdealGasModel::EulerIdealGasModel(const double _gamma)
:
    gamma(_gamma)
{}


/* Destructor */
EulerIdealGasModel::~EulerIdealGasModel()
{}

/* Return the ratio of specific heats */
const double EulerIdealGasModel::get_gamma() const
{
    return gamma;
}

/* Calculate the internal energy, e, from the
 * primitive variables */
double EulerIdealGasModel::prim_internal_energy(const Vec1D& prim) const
{
    // $p / (\rho (\gamma - 1))$
    return prim[3] / (prim[0] * (gamma - 1.0));
}

/* Calculate the total energy, E, from the
 * primitive variables */
double EulerIdealGasModel::prim_total_energy(const Vec1D& prim) const
{
    double e = prim_internal_energy(prim);
    
    // $\rho (\frac{1}{2} (u^2 + v^2) + e)$
    return prim[0] * (0.5 * (prim[1] * prim[1] + prim[2] * prim[2]) + e);
}

/* Calculate the enthalpy, h, from the
 * primitive variables */
double EulerIdealGasModel::prim_enthalpy(const Vec1D& prim) const
{
    double E = prim_total_energy(prim);
    
    // $(E + p) / \rho$
    return (E + prim[3]) / prim[0];
}

/* Calculate the local speed of sound, c, from the
 * primitive variables */
double EulerIdealGasModel::prim_speed_of_sound(const Vec1D& prim) const
{
    // $\sqrt(\gamma * p) / \rho$
    return sqrt(gamma * prim[3] / prim[0]);
}

/* Calculate the local Mach number from the
 * primitive variables */
double EulerIdealGasModel::prim_mach_number(const Vec1D& prim) const
{
    // $\sqrt(u^2 + v^2) / c$ 
    return sqrt(prim[1]*prim[1] + prim[2]*prim[2]) / prim_speed_of_sound(prim);
}
 
/* Convert the conservative flow variables to the
 * primitive flow variables */
Vec1D EulerIdealGasModel::cons_to_prim(const Vec1D& cons)
{
    double rho = cons[0];
    double rho_u = cons[1];
    double rho_v = cons[2];
    double E = cons[3];

    double u = rho_u / rho;
    double v = rho_v / rho;

    Vec1D prim(4);
    prim[0] = rho;
    prim[1] = u;
    prim[2] = v;
    prim[3] = (gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v));
    return prim;
}

/* Convert the primitive flow variables to the
 * conservative flow variables */
Vec1D EulerIdealGasModel::prim_to_cons(const Vec1D& prim)
{
    double rho = prim[0];
    double u = prim[1];
    double v = prim[2];
    double p = prim[3];

    Vec1D cons(4);
    cons[0] = rho;
    cons[1] = rho * u;
    cons[2] = rho * v;
    cons[3] = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
    return cons;
}

/* Calculate the conservative fluxes in the x-coordinate
 * direction from the conservative flow variables */
Vec1D EulerIdealGasModel::cons_flux_x(const Vec1D& cons)
{
    Vec1D prim = cons_to_prim(cons);
    double E = prim_total_energy(prim);

    double rho = prim[0];
    double u = prim[1];
    double v = prim[2];
    double p = prim[3];

    Vec1D flux(4);
    flux[0] = rho * u;
    flux[1] = rho * u * u + p;
    flux[2] = rho * u * v;
    flux[3] = u * (E + p);
    return flux;
}

/* Calculate the conservative fluxes in the y-coordinate
 * direction from the conservative flow variables */
Vec1D EulerIdealGasModel::cons_flux_y(const Vec1D& cons)
{
    Vec1D prim = cons_to_prim(cons);
    double E = prim_total_energy(prim);

    double rho = prim[0];
    double u = prim[1];
    double v = prim[2];
    double p = prim[3];

    Vec1D flux(4);
    flux[0] = rho * v;
    flux[1] = rho * u * v;
    flux[2] = rho * v * v + p;
    flux[3] = v * (E + p);
    return flux;
}

/* Rotate the flow values by an angle in degrees */
Vec1D EulerIdealGasModel::rotate_flow_values(const double angle_in_deg,
                                             const Vec1D& flow_values)
{
    double angle_in_rads = angle_in_deg * (M_PI/180.0);

    Vec1D rotated_values(flow_values.size());
    rotated_values[0] = flow_values[0];
    rotated_values[1] = cos(angle_in_rads) * flow_values[1] + sin(angle_in_rads) * flow_values[2];
    rotated_values[2] = -1.0 * sin(angle_in_rads) * flow_values[1] + cos(angle_in_rads) * flow_values[2];
    rotated_values[3] = flow_values[3];

    return rotated_values;
}

}

#endif
