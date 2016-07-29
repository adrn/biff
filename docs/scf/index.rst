
*********************************************
Self-consistent field expansions (`biff.scf`)
*********************************************

Introduction
============

For code blocks below, the following imports will be needed::

    >>> import astropy.units as u
    >>> import numpy as np
    >>> import biff.scf as bscf

Getting Started
===============

The two main uses of `biff.scf` are (1) to compute the expansion coefficients
given a continuous density distribution or discrete samples from a density
distribution, or (2) to evaluate the density, potential, and gradients of a
basis function expansion representation of a density distribution given a set of
coefficients.

TODO: compute coefficients

To evaluate properties of the potential or the density of the expansion
representation, Biff provides a :class:`~gala.potential.PotentialBase`
subclass--:class:`~biff.scf.SCFPotential`-- that supports all of the
functionality implemented in `gala` (especially `gala.potential`). To create an
instance, you must already have the expansion coefficients on hand (``Snlm`` for
the cosine terms and ``Tnlm`` for the sine terms)::

    >>>

.. toctree::
   :maxdepth: 1

   examples

API
===

.. automodapi:: biff.scf
