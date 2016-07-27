
*********************************************
Self-consistent field expansions (`biff.scf`)
*********************************************

Introduction
============

For code blocks below and any pages linked below, I assume the following
imports have already been excuted::

    >>> import astropy.units as u
    >>> import numpy as np
    >>> import biff.scf as bscf

Getting Started
===============

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
