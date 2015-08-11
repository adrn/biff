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
    double atan2(double y, double x) nogil
    double cos(double x) nogil
    double sin(double x) nogil

cdef extern from "gsl/gsl_sf_gamma.h":
    double gsl_sf_gamma(double x) nogil
    double gsl_sf_fact(unsigned int n) nogil

cdef extern from "src/coeff_helper.c":
    double RR_Plm_cosmphi(double r, double phi, double X, double a, int n, int l, int m) nogil

__all__ = ['Anlm_integrand']

cpdef Anlm_integrand(double phi, double X, double xsi,
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
