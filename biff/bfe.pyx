# coding: utf-8
# cython: boundscheck=False
# cython: nonecheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: profile=False
# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE=1

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

from collections import OrderedDict

from astropy.constants import G
import numpy as np
cimport numpy as np
from libc.math cimport M_PI

# Gala
from gala.units import galactic
from gala.potential.cpotential cimport CPotentialWrapper
from gala.potential.cpotential import CPotentialBase

cdef extern from "math.h":
    double sqrt(double x) nogil
    double atan2(double y, double x) nogil
    double cos(double x) nogil
    double sin(double x) nogil

cdef extern from "src/cpotential.h":
    enum:
        MAX_N_COMPONENTS = 16

    ctypedef double (*densityfunc)(double t, double *pars, double *q) nogil
    ctypedef double (*valuefunc)(double t, double *pars, double *q) nogil
    ctypedef void (*gradientfunc)(double t, double *pars, double *q, double *grad) nogil

    ctypedef struct CPotential:
        int n_components
        int n_dim
        densityfunc density[MAX_N_COMPONENTS]
        valuefunc value[MAX_N_COMPONENTS]
        gradientfunc gradient[MAX_N_COMPONENTS]
        int n_params[MAX_N_COMPONENTS]
        double *parameters[MAX_N_COMPONENTS]

cdef extern from "src/bfe_helper.h":
    double rho_nlm(double s, double phi, double X, int n, int l, int m) nogil
    double phi_nlm(double s, double phi, double X, int n, int l, int m) nogil
    double sph_grad_phi_nlm(double s, double phi, double X, int n, int l, int m, double *grad) nogil

cdef extern from "src/bfe.h":
    void scf_density_helper(double *xyz, int K, double M, double r_s,
                            double *Snlm, double *Tnlm,
                            int nmax, int lmax, double *dens) nogil
    void scf_potential_helper(double *xyz, int K, double G, double M, double r_s,
                              double *Snlm, double *Tnlm,
                              int nmax, int lmax, double *potv) nogil
    void scf_gradient_helper(double *xyz, int K, double G, double M, double r_s,
                             double *Snlm, double *Tnlm,
                             int nmax, int lmax, double *grad) nogil

    double scf_value(double t, double *pars, double *q) nogil
    double scf_density(double t, double *pars, double *q) nogil
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
        A 2D array of positions where ``axis=0`` are multiple positions
        and ``axis=1`` are the coordinate dimensions (x, y, z).
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

    scf_density_helper(&xyz[0,0], ncoords, M, r_s,
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
        A 2D array of positions where ``axis=0`` are multiple positions
        and ``axis=1`` are the coordinate dimensions (x, y, z).
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

    scf_potential_helper(&xyz[0,0], ncoords, G, M, r_s,
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
        A 2D array of positions where ``axis=0`` are multiple positions
        and ``axis=1`` are the coordinate dimensions (x, y, z).
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

    scf_gradient_helper(&xyz[0,0], ncoords, G, M, r_s,
                        &Snlm[0,0,0], &Tnlm[0,0,0],
                        nmax, lmax, &grad[0,0])

    return np.array(grad)

# ------------------------------------------------------

cdef class SCFWrapper(CPotentialWrapper):

    def __init__(self, G, parameters):
        cdef CPotential cp

        # This is the only code that needs to change per-potential
        cp.value[0] = <valuefunc>(scf_value)
        cp.density[0] = <densityfunc>(scf_density)
        cp.gradient[0] = <gradientfunc>(scf_gradient)
        # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        cp.n_components = 1
        self._params = np.array([G] + list(parameters), dtype=np.float64)
        self._n_params = np.array([len(self._params)], dtype=np.int32)
        cp.n_params = &(self._n_params[0])
        cp.parameters[0] = &(self._params[0])
        cp.n_dim = 3
        self.cpotential = cp

class SCFPotential(CPotentialBase):
    r"""
    SCFPotential(m, r_s, Snlm, Tnlm, units=None)

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
    def __init__(self, m, r_s, Snlm, Tnlm, units=None):
        Snlm = np.array(Snlm)
        Tnlm = np.array(Tnlm)

        parameters = OrderedDict()
        parameters['m'] = m
        parameters['r_s'] = r_s
        parameters['Snlm'] = Snlm
        parameters['Tnlm'] = Tnlm
        super(CPotentialBase, self).__init__(parameters, units=units)

        # specialized
        nmax = Tnlm.shape[0]-1
        lmax = Tnlm.shape[1]-1
        coeff = np.concatenate((Snlm.ravel(), Tnlm.ravel()))

        c_params = []
        c_params.append(self.parameters['m'].value)
        c_params.append(self.parameters['r_s'].value)
        c_params.append(nmax)
        c_params.append(lmax)
        c_params = c_params + coeff.tolist()
        self.c_parameters = np.array(c_params)
        self.c_instance = SCFWrapper(self.G, list(self.c_parameters))
