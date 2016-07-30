
*********************************************
Self-consistent field expansions (`biff.scf`)
*********************************************

Introduction
============

The two main uses of `biff.scf` are:

#. to compute the expansion coefficients given a continuous density distribution
   or discrete samples from a density distribution, or
#. to evaluate the density, potential, and gradients of a basis function
   expansion representation of a density distribution given a set of
   coefficients.

To compute expansion coefficients, the relevant functions are
`~biff.scf.compute_coeffs` and `~biff.scf.compute_coeffs_discrete`. This
implementation uses the notation from [L11]_: all expansion coefficients are
real, :math:`S_{nlm}` are the cosine coefficients, and :math:`T_{nlm}` are the
sine coefficients.

Once you have coefficients, there are two ways to evaluate properties of the
potential or the density of the expansion representation. Biff provides a
class-based interface :class:`~biff.scf.SCFPotential` that utilizes the
gravitational potential machinery implemented in `gala.potential` (and supports
all of the `gala` functionality, such as orbit integration and plotting). The
examples below use this interface.

As an alternate, there is also a functional interface to each relevant function:
`~biff.scf.density`, `~biff.scf.potential`, and `~biff.scf.gradient`.

Examples
--------
- :ref:`coeff-particle`
- :ref:`coeff-analytic`

.. toctree::
   :hidden:

   examples

API
===

.. automodapi:: biff.scf
