#ifndef __MODEL_H
#define __MODEL_H

#include "../../tensor/tensor.h"

namespace HyperFlow {

class Model
{
    
    public:
        
        /* Constructor */
        Model();

        /* Destructor */
        virtual ~Model();

        /* Calculate the local speed of sound from the
         * primitive variables */
        virtual double prim_speed_of_sound(const Vec1D& prim) const = 0;
 
        /* Convert the conservative flow variables to the
         * primitive flow variables */
        virtual Vec1D cons_to_prim(const Vec1D& cons) = 0;

        /* Convert the primitive flow variables to the
         * conservative flow variables */
        virtual Vec1D prim_to_cons(const Vec1D& prim) = 0;

        /* Calculate the conservative fluxes in the x-coordinate
         * direction from the conservative flow variables */
        virtual Vec1D cons_flux_x(const Vec1D& cons) = 0;

        /* Calculate the conservative fluxes in the y-coordinate
         * direction from the conservative flow variables */
        virtual Vec1D cons_flux_y(const Vec1D& cons) = 0;
};

}

#endif
