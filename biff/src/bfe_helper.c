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

void grad_phi_nlm(double r, double phi, double X, double r_s, int n, int l, int m, double *grad) {
    double dPhi_dr, dPhi_dphi, dPhi_dtheta;
    double s = r/r_s;

    // spherical coord stuff
    double sintheta = sqrt(1-X*X);
    double cosphi = cos(phi);
    double sinphi = sqrt(1-cosphi*cosphi);

    double cosmphi = cos(m*phi);
    double scoeff, leg, ggn;
    scoeff = pow(r,l) / pow(1+s,2*l+1) * gsl_sf_gegenpoly_n(n, 2*l + 1.5, (s-1)/(s+1));
    leg = gsl_sf_legendre_Plm(l, m, X);
    ggn = gsl_sf_gegenpoly_n(n, 1.5 + 2*l, (-1 + s)/(1 + s));

    // copied out of Mathematica
    if (n == 0) {
        dPhi_dr = -(pow(r_s,l)*pow(s,l-1)*pow(1+s,-2*l-3)*cosmphi*leg*
            (1 + s)*(l*(-1 + s) + s)*ggn);
    } else {
        dPhi_dr = -(pow(r_s,l)*pow(s,l-1)*pow(1+s,-2*l-3)*cosmphi*leg*
            (-2*(3 + 4*l)*s*gsl_sf_gegenpoly_n(n-1, 2.5 + 2*l, (-1 + s)/(1 + s)) +
            (1 + s)*(l*(-1 + s) + s)*ggn));
    }

    if (l==0) {
        dPhi_dtheta = 0.;
    } else if ((l-1) < m) {
        dPhi_dtheta = scoeff * leg * l*X*leg / (r*sintheta);
    } else {
        dPhi_dtheta = scoeff * leg * (-l*X*leg + (l+m)*gsl_sf_legendre_Plm(l-1, m, X)) / (-r*sintheta);
    }

    if (m == 0) {
        dPhi_dphi = 0.;
    } else {
        dPhi_dphi = scoeff * leg * (-m*sin(m*phi)) / (r*sintheta);
    }

    grad[0] = -(sintheta*cosphi*dPhi_dr + X*cosphi*dPhi_dtheta - sinphi*dPhi_dphi);
    grad[1] = -(sintheta*sinphi*dPhi_dr + X*sinphi*dPhi_dtheta + cosphi*dPhi_dphi);
    grad[2] = -(X*dPhi_dr - sintheta*dPhi_dtheta);
}
