#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_sf_gegenbauer.h"
#include <math.h>
#include "bfe_helper.h"

#define SQRT_FOURPI 3.544907701811031

double rho_nl(double s, int n, int l) {
    double RR, Knl;
    Knl = 0.5*n*(n+4*l+3) + (l+1)*(2*l+1);
    RR = Knl/(2*M_PI) * pow(s,l) / (s*pow(1+s,2*l+3)) * gsl_sf_gegenpoly_n(n, 2*l + 1.5, (s-1)/(s+1));
    return RR;
}
double rho_nlm(double s, double phi, double X, int n, int l, int m) {
    return rho_nl(s, n, l) * gsl_sf_legendre_Plm(l, m, X);
}

double phi_nl(double s, int n, int l) {
    return -pow(s,l) * pow(1+s, -2*l-1) * gsl_sf_gegenpoly_n(n, 2*l+1.5, (s-1)/(s+1));
}
double phi_nlm(double s, double phi, double X, int n, int l, int m) {
    return phi_nl(s, n, l) * gsl_sf_legendre_Plm(l, m, X);
}

void sph_grad_phi_nlm(double s, double phi, double X, int n, int l, int m, double *sphgrad) {
    double dPhi_dr, dPhi_dphi, dPhi_dtheta;

    // spherical coord stuff
    double sintheta = sqrt(1-X*X);
    double cosmphi = cos(m*phi);

    double Phi_nl, Ylm;
    Phi_nl = phi_nl(s, n, l);
    Ylm = gsl_sf_legendre_sphPlm(l, m, X);

    double ggn = gsl_sf_gegenpoly_n(n, 1.5 + 2*l, (-1 + s)/(1 + s));

    // copied out of Mathematica
    if (n == 0) {
        dPhi_dr = -(pow(s,l)*pow(s,l-1)*pow(1+s,-2*l-3)*cosmphi*Ylm*
            (1 + s)*(l*(-1 + s) + s)*ggn);
    } else {
        dPhi_dr = -(pow(s,l)*pow(s,l-1)*pow(1+s,-2*l-3)*cosmphi*Ylm*
            (-2*(3 + 4*l)*s*gsl_sf_gegenpoly_n(n-1, 2.5 + 2*l, (-1 + s)/(1 + s)) +
            (1 + s)*(l*(-1 + s) + s)*ggn));
    }

    if (l==0) {
        dPhi_dtheta = 0.;
    } else if ((l-1) < m) {
        dPhi_dtheta = l*X*Ylm / sintheta;
    } else {
        dPhi_dtheta = (l*X*Ylm - (l+m)*gsl_sf_legendre_sphPlm(l-1, m, X)) / sintheta;
    }
    dPhi_dtheta *= Phi_nl / s;

    if (m == 0) {
        dPhi_dphi = 0.;
    } else {
        dPhi_dphi = m / sintheta;
    }
    dPhi_dphi *= Ylm * Phi_nl / s;

    sphgrad[0] = SQRT_FOURPI * dPhi_dr;
    sphgrad[1] = dPhi_dtheta;
    sphgrad[2] = dPhi_dphi;
}
