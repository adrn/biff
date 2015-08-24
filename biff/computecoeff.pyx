# coding: utf-8
# cython: boundscheck=False
# cython: nonecheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: profile=False

""" THIS IS A THIN WRAPPER AROUND THE FUNCTIONS IN coeff_helper.c """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

import numpy as np
cimport numpy as np
from libc.math cimport M_PI

cdef extern from "math.h":
    double sqrt(double x) nogil
    double cos(double x) nogil
    double sin(double x) nogil

cdef extern from "src/coeff_helper.c":
    double c_Snlm_integrand(double phi, double X, double xsi, double density, int n, int l, int m)
    double c_Tnlm_integrand(double phi, double X, double xsi, double density, int n, int l, int m)

__all__ = ['Snlm_integrand', 'Tnlm_integrand']

cpdef Snlm_integrand(double phi, double X, double xsi,
                     density_func,
                     int n, int l, int m,
                     double M, double r_s):
    cdef:
        double s = (1 + xsi) / (1 - xsi)
        double r = s * r_s
        double x = r * cos(phi) * sqrt(1-X*X)
        double y = r * sin(phi) * sqrt(1-X*X)
        double z = r * X

    return c_Snlm_integrand(phi, X, xsi,
                            density_func(x, y, z, M, r_s) / M * r_s*r_s*r_s,
                            n, l, m)

cpdef Tnlm_integrand(double phi, double X, double xsi,
                     density_func,
                     int n, int l, int m,
                     double M, double r_s):
    cdef:
        double s = (1 + xsi) / (1 - xsi)
        double r = s * r_s
        double x = r * cos(phi) * sqrt(1-X*X)
        double y = r * sin(phi) * sqrt(1-X*X)
        double z = r * X

    return c_Tnlm_integrand(phi, X, xsi,
                            density_func(x, y, z, M, r_s) / M * r_s*r_s*r_s,
                            n, l, m)
