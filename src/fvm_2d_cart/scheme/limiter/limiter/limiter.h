#ifndef __LIMITER_H
#define __LIMITER_H

#include <algorithm>
#include <cmath>

namespace HyperFlow {

class Limiter
{

    public:

        /* Constructor with slope limiter parameters */
        Limiter(const double _omega = 0.0,
                const double _beta_im_half = 1.0,
                const double _beta_ip_half = 1.0,
                const double _tol = 1e-6);
    
        /* Destructor */
        virtual ~Limiter();
    
        /* Calculate the limit for a particular difference ratio */
        virtual double operator() (const double r) = 0;

    protected:
    
        /* Omega limiter parameter as per Toro, 2nd Ed. */
        double omega;
    
        /* Omega limiter parameter as per Toro, 2nd Ed. */
        double beta_im_half;
    
        /* Omega limiter parameter as per Toro, 2nd Ed.)*/
        double beta_ip_half;
    
        /* Tolerance for rounding up a small number to avoid 
         * division by zero */
        double tol;

        /* Ensure that the denominator is not set to zero
         * to avoid division-by-zero errors */
        virtual double denom_tol(const double denom);
    
        /* Calculate \xi_left as per Toro, 2nd Ed. */
        virtual double xi_left(const double r);
    
        /* Calculate \xi_right as per Toro, 2nd Ed. */
        virtual double xi_right(const double r);
    
        /* Calculate the minimum of \xi_left and \xi_right */
        virtual double min_xi(const double r);
};

}

#endif
