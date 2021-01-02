#ifndef __VANLEER_LIMITER_CPP
#define __VANLEER_LIMITER_CPP

#include "vanleer.h"

namespace HyperFlow {

/* Constructor */
VanLeerLimiter::VanLeerLimiter()
{}

/* Constructor with slope limiter parameters */
VanLeerLimiter::VanLeerLimiter(const double _omega,
                             const double _beta_im_half,
                             const double _beta_ip_half,
                             const double _tol)
:
    Limiter(_omega, 
            _beta_im_half,
            _beta_ip_half,
            _tol)
{}

/* Destructor */
VanLeerLimiter::~VanLeerLimiter()
{}

/* Calculate the limit for a particular difference ratio */
double VanLeerLimiter::operator()(const double r)
{
    if (r <= 0.0) {
        return 0.0;
    } else {
        return std::min(((2.0 * r) / (1.0 + r)), xi_right(r));
    }
}

}

#endif