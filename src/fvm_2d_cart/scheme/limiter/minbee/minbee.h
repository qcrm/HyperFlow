#ifndef __MINBEE_LIMITER_H
#define __MINBEE_LIMITER_H

#include "../limiter/limiter.h"

namespace HyperFlow {

class MinBeeLimiter
: 
	public Limiter
{

    public:
        
        /* Constructor */
        MinBeeLimiter();
        
        /* Constructor with slope limiter parameters */
        MinBeeLimiter(const double _omega = 0.0,
                      const double _beta_im_half = 1.0,
                      const double _beta_ip_half = 1.0,
                      const double _tol = 1e-6);
        
        /* Destructor */
        virtual ~MinBeeLimiter();
        
        /* Calculate the limit for a particular difference ratio */
        virtual double operator() (const double r);
};

}

#endif