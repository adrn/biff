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

double c_Anlm_integrand(double phi, double X, double xsi,
                        density_func,
                     int n, int l, int m,
                     double M, double r_s,
                     double[::1] args):
    """
    Anlm_integrand(phi, X, xsi, density_func, n, l, m, M, r_s, density_func_args)
    """
    cdef:
        double r = (1 + xsi) / (1 - xsi)
        double x = r * cos(phi) * sqrt(1-X*X)
        double y = r * sin(phi) * sqrt(1-X*X)
        double z = r * X

        double Knl, Inl, tmp, tmp2, krond

    Knl = 0.5*n*(n + 4*l + 3) + (l + 1)*(2*l + 1)
    tmp2 = (gsl_sf_gamma(n + 4*l + 3) / (gsl_sf_fact(n) *
            (n + 2*l + 1.5) * gsl_sf_gamma(2*l + 1.5)**2))
    if m == 0:
        krond = 1.
    else:
        krond = 0.

    Inl = (Knl / 2**(8*l+6) * tmp2 * (1 + krond) * M_PI *
           2/(2*l+1) * gsl_sf_fact(l+m) / gsl_sf_fact(l-m))
    tmp = 2. / ((1-xsi)*(1-xsi))
    return (RR_Plm_cosmphi(r, phi, X, r_s, n, l, m) * tmp / Inl *
            density_func(x, y, z, args) / M)
