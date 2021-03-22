#ifndef __EULER_IDEAL_GAS_H
#define __EULER_IDEAL_GAS_H

#include <cmath>
#include <memory>
#include <vector>

#include "../model/model.h"

namespace HyperFlow {

class EulerIdealGasModel : public Model
{
    
    public:

        /* Constructor */
        EulerIdealGasModel();

        /* Constructor via ratio of specific heats,
         * defaulting to air (\gamma = 1.4) */
        EulerIdealGasModel(const double _gamma = 1.4);

        /* Destructor */
        virtual ~EulerIdealGasModel();

        /* Return the ratio of specific heats */
        const double get_gamma() const;

        /* Calculate the internal energy, e, from the
         * primitive variables */
        double prim_internal_energy(const Vec1D& prim) const;
        
        /* Calculate the total energy from the
         * primitive variables */
        double prim_total_energy(const Vec1D& prim) const;
        
        /* Calculate the enthalpy from the
         * primitive variables */
        double prim_enthalpy(const Vec1D& prim) const;
        
        /* Calculate the local speed of sound from the
         * primitive variables */
        double prim_speed_of_sound(const Vec1D& prim) const;
        
        /* Calculate the local Mach number from the
         * primitive variables */
        double prim_mach_number(const Vec1D& prim) const;

        /* Convert the conservative flow variables to the
         * primitive flow variables */
        Vec1D cons_to_prim(const Vec1D& cons);

        /* Convert the primitive flow variables to the
         * conservative flow variables */
        Vec1D prim_to_cons(const Vec1D& prim);

        /* Calculate the conservative fluxes in the x-coordinate
         * direction from the conservative flow variables */
        Vec1D cons_flux_x(const Vec1D& cons);

        /* Calculate the conservative fluxes in the y-coordinate
         * direction from the conservative flow variables */
        Vec1D cons_flux_y(const Vec1D& cons);

        /* Rotate the flow values by an angle in degrees */
        Vec1D rotate_flow_values(const double angle_in_deg, const Vec1D& flow_values);
        
    private:

        /* Ratio of specific heats */
        double gamma;
};

}

#endif
