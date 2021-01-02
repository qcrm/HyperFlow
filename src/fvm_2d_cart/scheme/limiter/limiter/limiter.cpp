#ifndef __LIMITER_CPP
#define __LIMITER_CPP

#include "limiter.h"

namespace HyperFlow {

/* Constructor with slope limiter parameters */
Limiter::Limiter(const double _omega,
                 const double _beta_im_half,
                 const double _beta_ip_half,
                 const double _tol)
:
    omega(_omega),
    beta_im_half(_beta_im_half),
    beta_ip_half(_beta_ip_half),
    tol(_tol)
{}

/* Destructor */
Limiter::~Limiter()
{}

// Ensure that the denominator is not set to zero
// to avoid division-by-zero errors
double Limiter::denom_tol(const double denom)
{
    if (fabs(denom) < tol) {
        if (denom == 0.0) {
            return tol;
        }
        // Sign function
        return tol * ((denom > 0) - (denom < 0));
    } else {
        return denom;
    }
}

/* Calculate \xi_left as per Toro 2nd Ed. */
double Limiter::xi_left(const double r)
{
    double num = 2.0 * beta_im_half * r;
    double denom = 1.0 - omega + (1.0 + omega) * r;
    return num / denom_tol(denom);
}

/* Calculate \xi_right as per Toro 2nd Ed. */
double Limiter::xi_right(const double r)
{
    double num = 2.0 * beta_ip_half;
    double denom = 1.0 - omega + (1.0 + omega) * r;
    return num / denom_tol(denom);
}

/* Calculate the minimum of \xi_left and \xi_right */
double Limiter::min_xi(const double r)
{
    return std::min(xi_left(r), xi_right(r));
}

}

#endif