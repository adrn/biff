# coding: utf-8
# cython: boundscheck=False
# cython: nonecheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: profile=False

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

from astropy.constants import G
import numpy as np
cimport numpy as np
from libc.math cimport M_PI

# Gary
from gary.units import galactic
from gary.potential.cpotential cimport _CPotential
from gary.potential.cpotential import CPotentialBase

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

    double scf_value(double t, double *pars, double *q) nogil
    void scf_gradient(double t, double *pars, double *q, double *grad) nogil

__all__ = ['density', 'potential', 'gradient', 'SCFPotential']

cpdef density(double[:,::1] xyz,
              double[:,:,::1] Snlm, double[:,:,::1] Tnlm,
              int nmax, int lmax,
              double M=1., double r_s=1.):
    """
    density(xyz, Snlm, Tnlm, nmax, lmax, M=1, r_s=1)

    Compute the density of the basis function expansion
    at a set of positions given the expansion coefficients.

    Parameters
    ----------
    xyz : `~numpy.ndarray`
        A 2D array of positions where ``axis=0`` are the different
        points and ``axis=1`` are the coordinate dimensions x, y, z.
    Snlm : `~numpy.ndarray`
        A 3D array of expansion coefficients for the cosine terms
        of the expansion. This notation follows Lowing et al. (2011).
        The array should have shape ``(nmax+1,lmax+1,lmax+1)`` and any
        invalid terms (e.g., when m > l) will be ignored.
    Tnlm : `~numpy.ndarray`
        A 3D array of expansion coefficients for the sine terms
        of the expansion. This notation follows Lowing et al. (2011).
        The array should have shape ``(nmax+1,lmax+1,lmax+1)`` and any
        invalid terms (e.g., when m > l) will be ignored.
    nmax : int
        Number of radial expansion terms.
    lmax : int
        Maximum ``l`` value for the spherical harmonics.
    M : numeric (optional)
        Mass scale. Leave unset for dimensionless units.
    r_s : numeric (optional)
        Length scale. Leave unset for dimensionless units.

    Returns
    -------
    dens : `~numpy.ndarray`
        A 1D array of the density at each input position.
        Will have the same length as the input position array, ``xyz``.

    TODOs
    -----

    * Rework this so it accepts an array of nlm's with the same
        shape as the coefficients so the user doesn't always have to
        pass in a full 3D array.

    """

    cdef:
        int ncoords = xyz.shape[0]
        double[::1] dens = np.zeros(ncoords)

    c_density(&xyz[0,0], ncoords, M, r_s,
              &Snlm[0,0,0], &Tnlm[0,0,0],
              nmax, lmax, &dens[0])

    return np.array(dens)

cpdef potential(double[:,::1] xyz,
                double[:,:,::1] Snlm, double[:,:,::1] Tnlm,
                int nmax, int lmax,
                double G=1., double M=1., double r_s=1.):
    """
    potential(xyz, Snlm, Tnlm, nmax, lmax, G=1, M=1, r_s=1)

    Compute the gravitational potential of the basis function expansion
    at a set of positions given the expansion coefficients.

    Parameters
    ----------
    xyz : `~numpy.ndarray`
        A 2D array of positions where ``axis=0`` are the different
        points and ``axis=1`` are the coordinate dimensions x, y, z.
    Snlm : `~numpy.ndarray`
        A 3D array of expansion coefficients for the cosine terms
        of the expansion. This notation follows Lowing et al. (2011).
        The array should have shape ``(nmax+1,lmax+1,lmax+1)`` and any
        invalid terms (e.g., when m > l) will be ignored.
    Tnlm : `~numpy.ndarray`
        A 3D array of expansion coefficients for the sine terms
        of the expansion. This notation follows Lowing et al. (2011).
        The array should have shape ``(nmax+1,lmax+1,lmax+1)`` and any
        invalid terms (e.g., when m > l) will be ignored.
    nmax : int
        Number of radial expansion terms.
    lmax : int
        Maximum ``l`` value for the spherical harmonics.
    G : numeric (optional)
        Gravitational constant. Leave unset for dimensionless units.
    M : numeric (optional)
        Mass scale. Leave unset for dimensionless units.
    r_s : numeric (optional)
        Length scale. Leave unset for dimensionless units.

    Returns
    -------
    pot : `~numpy.ndarray`
        A 1D array of the value of the potential at each input position.
        Will have the same length as the input position array, ``xyz``.

    TODOs
    -----

    * Rework this so it accepts an array of nlm's with the same
        shape as the coefficients so the user doesn't always have to
        pass in a full 3D array.

    """
    cdef:
        int ncoords = xyz.shape[0]
        double[::1] potv = np.zeros(ncoords)

    c_potential(&xyz[0,0], ncoords, G, M, r_s,
                &Snlm[0,0,0], &Tnlm[0,0,0],
                nmax, lmax, &potv[0])

    return np.array(potv)

cpdef gradient(double[:,::1] xyz,
               double[:,:,::1] Snlm, double[:,:,::1] Tnlm,
               int nmax, int lmax,
               double G=1, double M=1, double r_s=1):
    """
    gradient(xyz, Snlm, Tnlm, nmax, lmax, G=1, M=1, r_s=1)

    Compute the gradient of the gravitational potential of the
    basis function expansion at a set of positions given the
    expansion coefficients.

    Parameters
    ----------
    xyz : `~numpy.ndarray`
        A 2D array of positions where ``axis=0`` are the different
        points and ``axis=1`` are the coordinate dimensions x, y, z.
    Snlm : `~numpy.ndarray`
        A 3D array of expansion coefficients for the cosine terms
        of the expansion. This notation follows Lowing et al. (2011).
        The array should have shape ``(nmax+1,lmax+1,lmax+1)`` and any
        invalid terms (e.g., when m > l) will be ignored.
    Tnlm : `~numpy.ndarray`
        A 3D array of expansion coefficients for the sine terms
        of the expansion. This notation follows Lowing et al. (2011).
        The array should have shape ``(nmax+1,lmax+1,lmax+1)`` and any
        invalid terms (e.g., when m > l) will be ignored.
    nmax : int
        Number of radial expansion terms.
    lmax : int
        Maximum ``l`` value for the spherical harmonics.
    G : numeric (optional)
        Gravitational constant. Leave unset for dimensionless units.
    M : numeric (optional)
        Mass scale. Leave unset for dimensionless units.
    r_s : numeric (optional)
        Length scale. Leave unset for dimensionless units.

    Returns
    -------
    grad : `~numpy.ndarray`
        A 2D array of the gradient of the potential at each input position.
        Will have the same shape as the input position array, ``xyz``.

    TODOs
    -----

    * Rework this so it accepts an array of nlm's with the same
        shape as the coefficients so the user doesn't always have to
        pass in a full 3D array.

    """
    cdef:
        int ncoords = xyz.shape[0]
        double[:,::1] grad = np.zeros((ncoords,3))

    c_gradient(&xyz[0,0], ncoords, G, M, r_s,
               &Snlm[0,0,0], &Tnlm[0,0,0],
               nmax, lmax, &grad[0,0])

    return np.array(grad)

# ------------------------------------------------------

cdef class _SCFPotential(_CPotential):
    # double[:,:,::1] sin_coeff, double[:,:,::1] cos_coeff):
    # np.ndarray[np.float64_t, ndim=3] sin_coeff,
    # np.ndarray[np.float64_t, ndim=3] cos_coeff):
    def __cinit__(self, double G, double m, double r_s,
                  int nmax, int lmax,
                  *args):
        self._parvec = np.concatenate([[G,m,r_s,nmax,lmax],args])
                                       # sin_coeff.ravel(),
                                       # cos_coeff.ravel()])
        self._parameters = &(self._parvec[0])
        self.c_value = &scf_value
        self.c_gradient = &scf_gradient

class SCFPotential(CPotentialBase):
    r"""
    SCFPotential(M, r_s, Snlm, Tnlm, units)

    An SCF / basis function expansion potential. Follows the
    convention used in Hernquist & Ostriker (1992) and
    Lowing et al. (2011) for representing all coefficients as
    real quantities.

    Parameters
    ----------
    m : numeric
        Scale mass.
    r_s : numeric
        Scale length.
    Snlm : array_like
        Array of coefficients for the cosine terms of the expansion.
        This should be a 3D array with shape `(nmax+1, lmax+1, lmax+1)`,
        where `nmax` is the number of radial expansion terms and `lmax`
        is the number of spherical harmonic `l` terms.
    Tnlm : array_like
        Array of coefficients for the sine terms of the expansion.
        This should be a 3D array with shape `(nmax+1, lmax+1, lmax+1)`,
        where `nmax` is the number of radial expansion terms and `lmax`
        is the number of spherical harmonic `l` terms.
    units : iterable
        Unique list of non-reducable units that specify (at minimum) the
        length, mass, time, and angle units.

    """
    def __init__(self, m, r_s, Snlm, Tnlm, units=galactic):
        Snlm = np.array(Snlm)
        Tnlm = np.array(Tnlm)
        self.G = G.decompose(units).value
        self.parameters = dict()
        self.parameters['m'] = m
        self.parameters['r_s'] = r_s
        self.parameters['Snlm'] = Snlm
        self.parameters['Tnlm'] = Tnlm
        super(SCFPotential, self).__init__(units=units)

        nmax = Tnlm.shape[0]-1
        lmax = Tnlm.shape[1]-1

        # c_params = self.parameters.copy()
        # c_params['G'] = self.G
        # c_params.pop('sin_coeff')
        # c_params.pop('cos_coeff')
        coeff = np.concatenate((Snlm.ravel(), Tnlm.ravel()))
        params1 = [self.G, self.parameters['m'], self.parameters['r_s'],
                   nmax, lmax]
        c_params = np.array(params1 + coeff.tolist())
        # self.c_instance = _SCFPotential(*coeff, **c_params)
        self.c_instance = _SCFPotential(*c_params)
