# coding: utf-8

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import numpy as np
import scipy.integrate as si

# Project
from ._computecoeff import Anlm_integrand

__all__ = ['compute_Anlm']

def compute_Anlm(density_func, nlm, M, r_s, args=(), **tplquad_kwargs):
    """
    Compute the expansion coefficients for representing the input
    density function as a basis function expansion.

    This is largely based on the formalism outlined in
    Hernquist & Ostriker (1992). The radial basis functions
    are taken to be Hernquist functions.

    Computing the coefficients involves computing triple integrals
    which are computationally expensive. For an example of how to
    parallelize the computation of the coefficients, see
    ``examples/parallel_compute_Anlm.py``.

    Parameters
    ----------
    density_func : function, callable
        A function or callable object to evaluate the density at a
        given position. The call format must be of the form:
        ``density_func(x, y, z, M, r_s, args)`` where ``x,y,z`` are
        cartesian coordinates, ``M`` is a scale mass, ``r_s`` a scale
        radius, and ``args`` is an iterable containing any other
        arguments needed by the density function.
    nlm : iterable
        A list or iterable of integers containing the values of n, l,
        and m, e.g., ``nlm = [n,l,m]``.
    M : numeric
        Scale mass.
    r_s : numeric
        Scale radius.
    args : iterable (optional)
        A list or iterable of any other arguments needed by the density
        function.

    tplquad_kwargs
        Any additional keyword arguments are passed through to
        `~scipy.integrate.tplquad`.

    Returns
    -------
    Anlm : float
        The value of the expansion coefficient.
    Anlm_err : float
        An estimate of the uncertainty in the coefficient value.

    """
    nlm = np.array(nlm).astype(np.int32)
    _args = np.array(args)

    Anlm = si.tplquad(Anlm_integrand,
                      -1., 1.,
                      lambda *args: -1., lambda *args: 1.,
                      lambda *args: 0., lambda *args: 2*np.pi,
                      args=(density_func, nlm[0], nlm[1], nlm[2], M, r_s, _args),
                      **tplquad_kwargs)

    return Anlm
