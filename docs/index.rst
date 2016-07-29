.. _biff:

****
Biff
****

Biff is a Python package for evaluating basis function expansions of mass
densities and gravitational potentials.

Currently, the only implemented expansion is the Self-Consistent Field (SCF)
method of Hernquist & Ostriker (1992; [HO92]_) using Hernquist radial functions
and spherical harmonics for angular functions. This implementation is based on
the formalism described in the original paper but using the notation of Lowing
et al. (2011; [L11]_).

.. note::

    Once it is more stable, Biff may be merged in to `gala.potential` .

Documentation
=============

.. toctree::
   :maxdepth: 1

   install
   scf

References
----------
.. [HO92] http://dx.doi.org/10.1086/171025
.. [L11] http://dx.doi.org/10.1111/j.1365-2966.2011.19222.x
