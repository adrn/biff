# coding: utf-8
# cython: boundscheck=False
# cython: nonecheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: profile=False

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

cdef extern from "src/bfe_helper.c":
    double phi_nlm(double r, double phi, double X, double r_s, int n, int l, int m) nogil
    double rho_nlm(double r, double phi, double X, double r_s, int n, int l, int m) nogil

__all__ = ['density']

cpdef density(double[:,::1] xyz,
              double G, double M, double r_s,
              double[:,:,::1] Anlm, int nmax, int lmax):
    """
    density()
    """

    cdef:
        int i,n,l,m
        int ncoords = xyz.shape[0]
        double r,X,phi
        double[::1] pot = np.zeros(ncoords)

    for i in range(ncoords):
        r = sqrt(xyz[i,0]*xyz[i,0] + xyz[i,1]*xyz[i,1] + xyz[i,2]*xyz[i,2])
        X = xyz[i,2]/r # cos(theta)
        phi = atan2(xyz[i,1], xyz[i,0])

        for n in range(nmax+1):
            for l in range(lmax+1):
                for m in range(l+1):
                    pot[i] += Anlm[n,l,m] * phi_nlm(r, phi, X, r_s, n, l, m)

        pot[i] *= G*M/r_s

    return np.array(pot)

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
