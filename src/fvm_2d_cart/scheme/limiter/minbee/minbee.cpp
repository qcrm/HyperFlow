#ifndef __MINBEE_LIMITER_CPP
#define __MINBEE_LIMITER_CPP

#include "minbee.h"

namespace HyperFlow {

/* Constructor */
MinBeeLimiter::MinBeeLimiter()
{}

/* Constructor with slope limiter parameters */
MinBeeLimiter::MinBeeLimiter(const double _omega,
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
MinBeeLimiter::~MinBeeLimiter()
{}

/* Calculate the limit for a particular difference ratio */
double MinBeeLimiter::operator()(const double r)
{
    if (r <= 0.0) {
        return 0.0;
    } else if ((0.0 < r) && (r <= 1.0)) {
        return r;
    } else {
        return std::min(1.0, xi_right(r));
    }
}

}

#endif