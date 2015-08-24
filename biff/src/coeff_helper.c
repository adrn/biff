#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_sf_gegenbauer.h"
#include "gsl/gsl_sf_gamma.h"
#include <math.h>
#include "coeff_helper.h"
#include "bfe_helper.h"
#include <complex.h>

#define SQRT_FOURPI 3.544907701811031

double cc_phi_nlm(double s, double phi, double X, int n, int l, int m) {
    /* Complex conjugate of phi_nlm */
    double tmp = pow(-1.,m) * gsl_sf_fact(l-m) / gsl_sf_fact(l+m);
    return phi_nl(s, phi, X, n, l) * tmp * gsl_sf_legendre_sphPlm(l, m, X);
}

double STnlm_integrand_help(double phi, double X, double xsi,
                            double density, int n, int l, int m) {
    /*
    Computes the integrand used to compute the expansion
    coefficients, Snlm, Tnlm. The integral is done over:

        * phi: azimuthal angle
        * X: cos(theta), where theta is the colatitude
            (e.g., from spherical coordinates typical to physicists)
        * xsi: (s-1)/(s+1), a radial coordinate mapped to the interval
            [-1,1] rather than [0,inf].
    */
    double _tmp = (1 - xsi);
    double s = (1 + xsi) / _tmp;

    // temporary variables
    double Knl, Inl, krond, _tmp2;

    Knl = 0.5*n*(n + 4*l + 3) + (l + 1)*(2*l + 1);
    if (m == 0) {
        krond = 1.;
    } else {
        krond = 0.;
    }

    _tmp2 = (gsl_sf_gamma(n + 4*l + 3) / (gsl_sf_fact(n) * (n + 2*l + 1.5) * pow(gsl_sf_gamma(2*l + 1.5),2)));
    Inl = (Knl / pow(2., 8*l+6) * _tmp2 * (1 + krond) * M_PI * 2./(2*l+1.) * gsl_sf_fact(l+m) / gsl_sf_fact(l-m));
    // Anl = (-pow(2., 8*l+6) / (Knl * 4*M_PI) *
    //        (gsl_sf_fact(n) * (n + 2*l + 1.5) * pow(gsl_sf_gamma(2*l + 1.5),2)) / gsl_sf_gamma(n+4*l+3));
    return -2. / (_tmp*_tmp) * s*s * cc_phi_nlm(s, phi, X, n, l, m) / Inl * density;
}

extern double c_Snlm_integrand(double phi, double X, double xsi,
                               double density, int n, int l, int m) {
    return STnlm_integrand_help(phi, X, xsi, density, n, l, m) * cos(m*phi);
}

extern double c_Tnlm_integrand(double phi, double X, double xsi,
                               double density, int n, int l, int m) {
    return STnlm_integrand_help(phi, X, xsi, density, n, l, m) * sin(m*phi);
}
