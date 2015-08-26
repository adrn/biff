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
    double rho_nlm(double s, double phi, double X, int n, int l, int m) nogil
    double phi_nlm(double s, double phi, double X, int n, int l, int m) nogil
    double sph_grad_phi_nlm(double s, double phi, double X, int n, int l, int m, double *grad) nogil

cdef extern from "src/bfe.c":
    void c_density(double *xyz, int K, double M, double r_s,
                   double *Snlm, double *Tnlm,
                   int nmax, int lmax, double *dens) nogil
    void c_potential(double *xyz, int K, double G, double M, double r_s,
                     double *Snlm, double *Tnlm,
                     int nmax, int lmax, double *potv) nogil
    void c_gradient(double *xyz, int K, double G, double M, double r_s,
                    double *Snlm, double *Tnlm,
                    int nmax, int lmax, double *grad) nogil

__all__ = ['density', 'potential', 'gradient']

cpdef density(double[:,::1] xyz,
              double M, double r_s,
              double[:,:,::1] Snlm, double[:,:,::1] Tnlm,
              int nmax, int lmax):

    cdef:
        int ncoords = xyz.shape[0]
        double[::1] dens = np.zeros(ncoords)

    c_density(&xyz[0,0], ncoords, M, r_s,
              &Snlm[0,0,0], &Tnlm[0,0,0],
              nmax, lmax, &dens[0])

    return np.array(dens)

cpdef potential(double[:,::1] xyz,
                double G, double M, double r_s,
                double[:,:,::1] Snlm, double[:,:,::1] Tnlm,
                int nmax, int lmax):

    cdef:
        int ncoords = xyz.shape[0]
        double[::1] potv = np.zeros(ncoords)

    c_potential(&xyz[0,0], ncoords, G, M, r_s,
                &Snlm[0,0,0], &Tnlm[0,0,0],
                nmax, lmax, &potv[0])

    return np.array(potv)

cpdef gradient(double[:,::1] xyz,
               double G, double M, double r_s,
               double[:,:,::1] Snlm, double[:,:,::1] Tnlm,
               int nmax, int lmax):

    cdef:
        int ncoords = xyz.shape[0]
        double[:,::1] grad = np.zeros((ncoords,3))

    c_gradient(&xyz[0,0], ncoords, G, M, r_s,
               &Snlm[0,0,0], &Tnlm[0,0,0],
               nmax, lmax, &grad[0,0])

    return np.array(grad)
