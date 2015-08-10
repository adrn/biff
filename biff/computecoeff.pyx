# coding: utf-8
# cython: boundscheck=False
# cython: nonecheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: profile=False

""" Basis Function Expansion in Cython """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

import numpy as np
cimport numpy as np
from libc.math cimport M_PI

cdef extern from "math.h":
    double sqrt(double x) nogil

cdef extern from "src/coeff_helper.c":
    ctypedef double (*DensityFunc)(double, double, double, double, double, double*) nogil
    double Anlm_integrand(DensityFunc density, double phi, double X, double xsi,
                          int n, int l, int m, double M, double r_s, double *args) nogil

__all__ = ['compute_Anlm']

cdef double test_density(double x, double y, double z,
                         double M, double *args) nogil:
    cdef double r = sqrt(x*x + y*y + z*z)
    return M/(2.*M_PI) * args[0] / (r*(r + args[0])**3)

cpdef compute_Anlm(density_func):
    cdef:
        double M = 1E10
        double c = 2.
        double[::1] args = np.array([M, c])

    Anlm_integrand(<DensityFunc>density_func,
                   0., 0., 0.,
                   0, 0, 0,
                   M, c,
                   &args[0])














# cpdef phi_nlm(double r, double phi, double costheta, int n, int l, int m):
#     cdef double A,B
#     A = s**l / (1+s)**(2*l+1) * gsl_sf_gegenpoly_n(n, 2*l + 1.5, (s-1)/(s+1))
#     B = gsl_sf_legendre_Plm(l, m, costheta)
#     return -A * B * cos(m*phi)

# cpdef rho_nlm(double s, double phi, double costheta, int n, int l, int m):
#     cdef double A,B,Knl
#     Knl = 0.5*n*(n+4*l+3) + (l+1)*(2*l+1)
#     A = Knl/(2*M_PI) * s**l / (s*(1+s)**(2*l+3)) * gsl_sf_gegenpoly_n(n, 2*l + 1.5, (s-1)/(s+1))
#     B = gsl_sf_legendre_Plm(l, m, costheta)
#     return A * B * cos(m*phi)

# cpdef apw_value(double[::1] xyz,
#                 double G, double M, double r_s,
#                 double[:,:,::1] Anlm, int nmax, int lmax):
#     cdef:
#         double r = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2])
#         double s = r / r_s
#         double costheta = xyz[2]/r
#         double phi = atan2(xyz[1], xyz[0])
#         int n,l,m
#         double pot = 0.

#     for n in range(nmax+1):
#         for l in range(lmax+1):
#             for m in range(-l,l+1):
#                 pot += Anlm[n,l,m] * phi_nlm(s, phi, costheta, n, l, m)

#     return G*M/r_s * pot

# cpdef apw_density(double[::1] xyz,
#                   double G, double M, double r_s,
#                   double[:,:,::1] Anlm, int nmax, int lmax):
#     cdef:
#         double r = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2])
#         double s = r / r_s
#         double costheta = xyz[2]/r
#         double phi = atan2(xyz[1], xyz[0])
#         int n,l,m
#         double dens = 0.

#     for n in range(nmax+1):
#         for l in range(lmax+1):
#             for m in range(-l,l+1):
#                 dens += Anlm[n,l,m] * rho_nlm(s, phi, costheta, n, l, m)

#     return M/(4*M_PI*r_s*r_s*r_s) * dens
