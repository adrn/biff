
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

Example: Computing expansion coefficients from particle positions
-----------------------------------------------------------------

TODO: Samples from a density distribution...

Example: Computing expansion coefficients for an analytic density
-----------------------------------------------------------------

Here's another example of integrating the
`Lorenz equations <https://en.wikipedia.org/wiki/Lorenz_system>`_, a 3D
nonlinear system::

    >>> def F(t,w,sigma,rho,beta):
    ...     x,y,z,px,py,pz = w
    ...     return np.array([sigma*(y-x), x*(rho-z)-y, x*y-beta*z, 0., 0., 0.]).reshape(w.shape)
    >>> sigma, rho, beta = 10., 28., 8/3.
    >>> integrator = gi.DOPRI853Integrator(F, func_args=(sigma, rho, beta))
    >>> orbit = integrator.run([0.5,0.5,0.5,0,0,0], dt=1E-2, n_steps=1E4)
    >>> fig = orbit.plot()

.. plot::
    :align: center

    import astropy.units as u
    import matplotlib.pyplot as pl
    import numpy as np
    import gala.integrate as gi

API
===

.. automodapi:: biff.scf
