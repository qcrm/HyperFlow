#ifndef __HLLC_EULER_RIEMANN_SOLVER_CPP
#define __HLLC_EULER_RIEMANN_SOLVER_CPP

#include "riemann_hllc.h"

namespace HyperFlow {

/* Default Constructor */
HLLCEulerRiemannSolver::HLLCEulerRiemannSolver()
{};

/* Constructor from the Euler model */
HLLCEulerRiemannSolver::HLLCEulerRiemannSolver(
    const std::shared_ptr<EulerIdealGasModel>& _model
)
:
    model(_model)
{}

/* Destructor */
HLLCEulerRiemannSolver::~HLLCEulerRiemannSolver()
{}

/* Evaluate Riemann problem */
Vec1D HLLCEulerRiemannSolver::operator() (const Vec1D& cons_left,
                                          const Vec1D& cons_right,
                                          const Direction& dir)
{
    Vec1D hllc_flux = calc_hllc_flux(cons_left, cons_right, dir);
    return hllc_flux;
}

/* Calculate the velocities in the starred region on
 * either side of the contact discontinuity,
 * U_L_star and U_R_star
 * See Toro, 2nd Ed., pg 323, eqn 10.33. */
Vec1D HLLCEulerRiemannSolver::calc_U_k(const Vec1D& wK, 
                                       const double S_K,
                                       const double S_star,
                                       const Direction& dir) {
    double dK = wK[0];
    double uK = wK[1];
    double vK = wK[2];
    double pK = wK[3];
    double EK = model->prim_total_energy(wK);

    Vec1D U_K(4);

    if (dir == Direction::x) {
        double coef = dK * ((S_K - uK) / (S_K - S_star));
        U_K[0] = coef;
        U_K[1] = coef * S_star;
        U_K[2] = coef * vK;
        U_K[3] = coef * (EK/dK + (S_star - uK) * (S_star + (pK / (dK * (S_K - uK)))));
    } else {
        double coef = dK * ((S_K - vK) / (S_K - S_star));
        U_K[0] = coef;
        U_K[1] = coef * uK;
        U_K[2] = coef * S_star;
        U_K[3] = coef * (EK/dK + (S_star - vK) * (S_star + (pK / (dK * (S_K - vK)))));
    }
    return U_K;
}

/* Calculate the wave speed S_{*} in the starred region.
 * See Toro, 2nd Ed., pg 327, eqn 10.58. */
double HLLCEulerRiemannSolver::calc_S_star(const Vec1D& wL,
                                           const Vec1D& wR,
                                           const double S_L,
                                           const double S_R,
                                           const Direction& dir) {
    double dL = wL[0];
    double uL = wL[1];
    double vL = wL[2];
    double pL = wL[3];

    double dR = wR[0];
    double uR = wR[1];
    double vR = wR[2];
    double pR = wR[3];

    double velL = 0.0;
    double velR = 0.0;

    if (dir == Direction::x) {
        velL = uL;
        velR = uR;
    } else {
        velL = vL;
        velR = vR;
    }

    double num = pR - pL + dL * velL * (S_L - velL) - dR * velR * (S_R - velR);
    double denom = dL * (S_L - velL) - dR * (S_R - velR);
    
    return num/denom;
}

/* Calculate the wave speed estimates directly.
 * See Toro, 2nd Ed., pg 324. */
Vec1D HLLCEulerRiemannSolver::calc_wave_speeds_direct(const Vec1D& wL,
                                                      const Vec1D& wR,
                                                      const Direction& dir) {
    Vec1D waves(2, 0.0);
    
    double velL = 0.0;
    double velR = 0.0;

    if (dir == Direction::x) {
        velL = wL[1];
        velR = wR[1];
    } else {
        velL = wL[2];
        velR = wR[2];
    }

    waves[0] = velL - sqrt(model->get_gamma() * wL[3] / wL[0]);
    waves[1] = velR + sqrt(model->get_gamma() * wR[3] / wR[0]);

    return waves;
}

/* Calculate the wave speed estimates directly using
 * an alternative estimation mechanism. This improves
 * oscillations in the solution.
 * See Toro, 2nd Ed., pg 324, eqn 10.38. */
Vec1D HLLCEulerRiemannSolver::calc_wave_speeds_direct_min_max(const Vec1D& wL,
                                                              const Vec1D& wR,
                                                              const Direction& dir) {
    double a_L = sqrt(model->get_gamma() * wL[3] / wL[0]);
    double a_R = sqrt(model->get_gamma() * wR[3] / wR[0]);

    Vec1D waves(2, 0.0);

    double velL = 0.0;
    double velR = 0.0;

    if (dir == Direction::x) {
        velL = wL[1];
        velR = wR[1];
    } else {
        velL = wL[2];
        velR = wR[2];
    }

    waves[0] = std::min(velL - a_L, velR - a_R);
    waves[1] = std::max(velL + a_L, velR + a_R);
    
    return waves;
}

/* Calculate the fluxes in the starred region. */
Vec1D HLLCEulerRiemannSolver::calc_F_star(const Vec1D& FK,
                                          const Vec1D& U_K_star,
                                          const Vec1D& U_K,
                                          const double S_K) {
    Vec1D F_K_star(4);
    
    for (unsigned int i=0; i<4; i++) {
        F_K_star[i] = FK[i] + S_K * (U_K_star[i] - U_K[i]);
    }
    
    return F_K_star;
}

/* Calculate the Harten-Lax-van Leer-Contact (HLLC) of
 * Toro to obtain the Euler fluxes.
 * See Toro, 2nd Ed., pg 323. */
Vec1D HLLCEulerRiemannSolver::calc_hllc_flux(const Vec1D& U_L,
                                             const Vec1D& U_R,
                                             const Direction& dir) {
    Vec1D wL = model->cons_to_prim(U_L);
    Vec1D wR = model->cons_to_prim(U_R);

    Vec1D waves = calc_wave_speeds_direct_min_max(wL, wR, dir);

    double S_L = waves[0];
    double S_R = waves[1];

    Vec1D F_L;
    Vec1D F_R;

    if (dir == Direction::x) {
        F_L = model->cons_flux_x(U_L);
        F_R = model->cons_flux_x(U_R);
    } else {
        F_L = model->cons_flux_y(U_L);
        F_R = model->cons_flux_y(U_R);
    }

    double S_star = calc_S_star(wL, wR, S_L, S_R, dir);
    
    Vec1D U_L_star = calc_U_k(wL, S_L, S_star, dir);
    Vec1D U_R_star = calc_U_k(wR, S_R, S_star, dir);

    Vec1D F_L_star = calc_F_star(F_L, U_L_star, U_L, S_L);
    Vec1D F_R_star = calc_F_star(F_R, U_R_star, U_R, S_R);

    if (0.0 <= S_L) {
        return F_L;
    } else if (S_L <= 0.0 && 0.0 <= S_star) {
        return F_L_star;
    } else if (S_star <= 0.0 && 0.0 <= S_R) {
        return F_R_star;
    } else if (0.0 >= S_R) {
        return F_R;
    } else {
        /* Helpful debug information for inconsistent states */
        std::cout << "HLLC: ERROR - Should not get here." << std::endl;
        return Vec1D(4, 0.0);
    }
}

}

#endif
