***************
Biff (``biff``)
***************

Biff is a Python package for evaluating basis function expansions of mass
densities and gravitational potentials using the Self-Consistent Field (SCF)
method of Hernquist & Ostriker (1992) [1]_.

This is largely based on the formalism described in Hernquist & Ostriker (1992)
but using the notation of Lowing et al. (2011). The radial basis functions are
taken to be Hernquist functions and the angular functions are spherical
harmonics.

Installation
============


Overview
========

TODO: note difference in coefficients from HO92 - use Lowing form

TODO: Common thing to do for full grid is:

    >>> nmax = 6
    >>> lmax = 4
    >>> nlm = np.mgrid[0:nmax+1, 0:lmax+1, 0:lmax+1]
    >>> nlm.shape
    (3, 7, 5, 5)
    >>> nlm.reshape(3,-1).T.shape
    (175, 3)


Reference/API
=============

.. automodapi:: biff
