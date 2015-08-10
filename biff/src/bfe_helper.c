#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_sf_gegenbauer.h"
#include <math.h>
#include "bfe_helper.h"

double phi_nlm(double r, double phi, double X, double r_s, int n, int l, int m) {
    double A,B;
    double s = r/r_s;
    A = pow(r,l) / pow(1+s,2*l+1) * gsl_sf_gegenpoly_n(n, 2*l + 1.5, (s-1)/(s+1));
    B = gsl_sf_legendre_Plm(l, m, X);
    return -A * B * cos(m*phi);
}

double rho_nlm(double r, double phi, double X, double r_s, int n, int l, int m) {
    double A,B,Knl;
    double s = r/r_s;
    Knl = 0.5*n*(n+4*l+3) + (l+1)*(2*l+1);
    A = Knl/(2*M_PI) * pow(r,l) / (s*pow(1+s,2*l+3)) * gsl_sf_gegenpoly_n(n, 2*l + 1.5, (s-1)/(s+1));
    B = gsl_sf_legendre_Plm(l, m, X);
    return A * B * cos(m*phi);
}
