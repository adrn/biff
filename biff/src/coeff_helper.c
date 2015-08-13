#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_sf_gegenbauer.h"
#include <math.h>
#include "coeff_helper.h"
#include <complex.h>

double cc_phi_nlm(double s, double phi, double X, int n, int l, int m) {
    /* Complex conjugate of phi_nlm */
    double RR;
    double sqrt_fourpi = 3.544907701811031;
    RR = -pow(s,l) * pow(1+s, -2*l-1) * gsl_sf_gegenpoly_n(n, 2*l+1.5, (s-1)/(s+1));
    return sqrt_fourpi * RR * gsl_sf_legendre_sphPlm(l, m, X);
}

double STnlm_integrand_help(double phi, double X, double xsi,
                            double density, int n, int l, int m) {
    /*
    Computes the integrand used to compute the expansion
    coefficients, Anlm. The integral is done over:

        * phi: azimuthal angle
        * X: cos(theta), where theta is the colatitude
            (e.g., from spherical coordinates typical to physicists)
        * xsi: (s-1)/(s+1), a radial coordinate mapped to the interval
            [-1,1] rather than [0,inf].
    */
    double _tmp = (1 - xsi);
    double s = (1 + xsi) / _tmp;

    // temporary variables
    double Knl, Anl, krond;

    Knl = 0.5*n*(n + 4*l + 3) + (l + 1)*(2*l + 1);
    if (m == 0) {
        krond = 1.;
    } else {
        krond = 0.;
    }

    Anl = (-pow(2., 8*l+6) / (Knl * 4*M_PI) *
           (gsl_sf_fact(n) * (n + 2*l + 1.5) * pow(gsl_sf_gamma(2*l + 1.5),2)) / gsl_sf_gamma(n+4*l+3));
    return Anl * 2 / (_tmp*_tmp) * s*s * cc_phi_nlm(s, phi, X, n, l, m);
}

extern double Snlm_integrand(double phi, double X, double xsi,
                             double density, int n, int l, int m) {
    return STnlm_integrand_help(phi, X, xsi, density, n, l, m) * cos(m*phi);
}

extern double Tnlm_integrand(double phi, double X, double xsi,
                             double density, int n, int l, int m) {
    return STnlm_integrand_help(phi, X, xsi, density, n, l, m) * sin(m*phi);
}
