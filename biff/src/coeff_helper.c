#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_sf_gegenbauer.h"
#include <math.h>
#include "coeff_helper.h"

double RR_Plm_cosmphi(double r, double phi, double costheta,
                      double r_s, int n, int l, int m) {
    /*
    Computes a portion of the integrand for computing the
    expansion coefficients, Anlm.
    */
    double RR;
    double s = r/r_s;
    RR = pow(r,l) * pow(1+s,-2*l-1) * gsl_sf_gegenpoly_n(n, 2*l + 1.5, (s-1)/(s+1));
    return RR * gsl_sf_legendre_Plm(l, m, costheta) * cos(m*phi) * r*r;
}
