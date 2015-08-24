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
    double phi_nlm(double s, double phi, double X, int n, int l, int m) nogil
    double rho_nlm(double s, double phi, double X, int n, int l, int m) nogil
    double sph_grad_phi_nlm(double s, double phi, double X, int n, int l, int m, double *grad) nogil

__all__ = ['density', 'potential', 'gradient']

cpdef density(double[:,::1] xyz,
              double M, double r_s,
              double[:,:,::1] Snlm, double[:,:,::1] Tnlm,
              int nmax, int lmax):
    """
    density(xyz, M, r_s, Snlm, Tnlm, nmax, lmax)
    """

    cdef:
        int i,n,l,m
        int ncoords = xyz.shape[0]
        double r,s,X,phi
        double[::1] dens = np.zeros(ncoords)

    for i in range(ncoords):
        r = sqrt(xyz[i,0]*xyz[i,0] + xyz[i,1]*xyz[i,1] + xyz[i,2]*xyz[i,2])
        s = r/r_s
        X = xyz[i,2]/r # cos(theta)
        phi = atan2(xyz[i,1], xyz[i,0])

        for n in range(nmax+1):
            for l in range(lmax+1):
                for m in range(l+1):
                    dens[i] += rho_nlm(s, phi, X, n, l, m) * (Snlm[n,l,m]*cos(m*phi) +
                                                              Tnlm[n,l,m]*sin(m*phi))

        dens[i] *= M/(r_s*r_s*r_s)

    return np.array(dens)

cpdef potential(double[:,::1] xyz,
                double G, double M, double r_s,
                double[:,:,::1] Snlm, double[:,:,::1] Tnlm,
                int nmax, int lmax):
    """
    potential(xyz, G, M, r_s, Snlm, Tnlm, nmax, lmax)
    """

    cdef:
        int i,n,l,m
        int ncoords = xyz.shape[0]
        double r,s,X,phi
        double[::1] pot = np.zeros(ncoords)

    for i in range(ncoords):
        r = sqrt(xyz[i,0]*xyz[i,0] + xyz[i,1]*xyz[i,1] + xyz[i,2]*xyz[i,2])
        s = r/r_s
        X = xyz[i,2]/r # cos(theta)
        phi = atan2(xyz[i,1], xyz[i,0])

        for n in range(nmax+1):
            for l in range(lmax+1):
                for m in range(l+1):
                    pot[i] += rho_nlm(s, phi, X, n, l, m) * (Snlm[n,l,m]*cos(m*phi) +
                                                             Tnlm[n,l,m]*sin(m*phi))

        pot[i] *= G*M/r_s

    return np.array(pot)

cpdef gradient(double[:,::1] xyz,
               double G, double M, double r_s,
               double[:,:,::1] Snlm, double[:,:,::1] Tnlm,
               int nmax, int lmax):
    """
    gradient(xyz, G, M, r_s, Snlm, Tnlm, nmax, lmax)
    """

    cdef:
        int i,n,l,m
        int ncoords = xyz.shape[0]
        double r,s,X,phi
        double[:,::1] grad = np.zeros((ncoords,3))
        double[::1] tmp_grad = np.zeros(3)
        double tmp, tmp2, sintheta, cosphi, sinphi

    for i in range(ncoords):
        r = sqrt(xyz[i,0]*xyz[i,0] + xyz[i,1]*xyz[i,1] + xyz[i,2]*xyz[i,2])
        s = r/r_s
        X = xyz[i,2]/r # cos(theta)
        phi = atan2(xyz[i,1], xyz[i,0])

        sintheta = sqrt(1 - X*X)
        cosphi = cos(phi)
        sinphi = sin(phi)

        for n in range(nmax+1):
            for l in range(lmax+1):
                for m in range(l+1):
                    tmp = (Snlm[n,l,m]*cos(m*phi) + Tnlm[n,l,m]*sin(m*phi))
                    tmp2 = (Tnlm[n,l,m]*cos(m*phi) - Snlm[n,l,m]*sin(m*phi))

                    sph_grad_phi_nlm(s, phi, X, n, l, m, &tmp_grad[0])
                    grad[i,0] += (sintheta*cosphi*tmp_grad[0] + X*cosphi*tmp_grad[1] - sinphi*tmp_grad[2]) * tmp
                    grad[i,1] += (sintheta*sinphi*tmp_grad[0] + X*sinphi*tmp_grad[1] + cosphi*tmp_grad[2]) * tmp
                    grad[i,2] += (X*tmp_grad[0] - sintheta*tmp_grad[1]) * tmp2

        grad[i,0] *= G*M/(r_s*r_s)
        grad[i,1] *= G*M/(r_s*r_s)
        grad[i,2] *= G*M/(r_s*r_s)

    return np.array(grad)
