#ifndef __HLLC_EULER_RIEMANN_SOLVER_H
#define __HLLC_EULER_RIEMANN_SOLVER_H

#include <algorithm>
#include <cmath>
#include <iostream>

#include "../riemann/riemann.h"
#include "../../../model/euler_ideal_gas/euler_ideal_gas.h"

namespace HyperFlow {

class HLLCEulerRiemannSolver
:
    public RiemannSolver
{

    public:
    
        /* Constructor */
        HLLCEulerRiemannSolver(); 
        
        /* Construct from the Euler model */
        HLLCEulerRiemannSolver(const std::shared_ptr<EulerIdealGasModel>& _model);
        
        /* Destructor */
        virtual ~HLLCEulerRiemannSolver();
        
        /* Evaluate Riemann problem */ 
        virtual Vec1D operator() (const Vec1D& cons_left,
                                  const Vec1D& cons_right,
                                  const Direction& dir);

    private:
        
        /* Calculate the velocities in the starred region on
         * either side of the contact discontinuity,
         * U_L_star and U_R_star
         * See Toro, 2nd Ed., pg 323, eqn 10.33. */
        Vec1D calc_U_k(const Vec1D& wK, 
                       const double S_K,
                       const double S_star,
                       const Direction& dir);
        
        /* Calculate the wave speed S_{*} in the starred region.
         * See Toro, 2nd Ed., pg 327, eqn 10.58. */
        double calc_S_star(const Vec1D& wL,
                           const Vec1D& wR,
                           const double S_L,
                           const double S_R,
                           const Direction& dir);
        
        /* Calculate the wave speed estimates directly.
         * See Toro, 2nd Ed., pg 324. */
        Vec1D calc_wave_speeds_direct(const Vec1D& wL,
                                      const Vec1D& wR,
                                      const Direction& dir);
        
        /* Calculate the wave speed estimates directly using
         * an alternative estimation mechanism. This improves
         * oscillations in the solution.
         * See Toro, 2nd Ed., pg 324, eqn 10.38. */
        Vec1D calc_wave_speeds_direct_min_max(const Vec1D& wL,
                                              const Vec1D& wR,
                                              const Direction& dir);
        
        /* Calculate the fluxes in the starred region. */
        Vec1D calc_F_star(const Vec1D& FK,
                          const Vec1D& U_K_star,
                          const Vec1D& U_K,
                          const double S_K);
        
        /* Calculate the Harten-Lax-van Leer-Contact (HLLC) of
         * Toro to obtain the Euler fluxes.
         * See Toro, 2nd Ed., pg 323. */
        Vec1D calc_hllc_flux(const Vec1D& U_L,
                             const Vec1D& U_R,
                             const Direction& dir);

        /* Pointer to the model equations */
        std::shared_ptr<EulerIdealGasModel> model;
        
};

}

#endif
