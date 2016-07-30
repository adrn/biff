********
Examples
********

For code blocks below, the following imports will be needed::

    >>> import astropy.units as u
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import biff.scf as bscf

.. _coeff-particle:

Computing expansion coefficients from particle positions
--------------------------------------------------------

To compute expansion coefficients for a distribution of particles or discrete
samples from a density distribution, use `~biff.scf.compute_coeffs_discrete`. In
this example, we will generate particle positions from a Plummer density
profile, compute the expansion coefficients assuming spherical symmetry, then
re-compute the expanion coefficients and variances (Weinberg 1996; [W96]_)
allowing for non-spherical terms (e.g., :math:`l,m>0`).

We'll start by generating samples from a Plummer sphere (see Section 3 of
[HMV11]_ for more details). To do this, we will use inverse transform sampling
by inverting the cumulative mass function (in this case, the mass enclosed):

.. math::

    \rho(r) &= \frac{M}{\frac{4}{3}\pi a^3} \, \left(1 + \frac{r^2}{a^2}\right)^{-5/2}

    m(<r) &= \frac{M \, r^3}{(r^2 + a^2)^{3/2}}

    r(\mu) &= a \, (\mu^{-2/3} - 1)^{-1/2}

    \mu &= m(<r) / M

For simplicity, we will work with units in which :math:`a=1` and :math:`M=1`. To
generate radii, we first randomly generate values of :math:`\mu` uniformly
distributed between 0 and 1, then compute the value of :math:`r` for each
sample; the radii will then be distributed following a Plummer profile. For this
example, we'll use 16384 samples::

    def sample_r(size=1):
        mu = np.random.random(size=size)
        return 1 / np.sqrt(mu**(-2/3) - 1)

    n_samples = 16384
    r = sample_r(size=n_samples)

Let's plot the density profile derived from these samples vs. the true profile:

.. plot::
    :align: center
    :context:

    from astropy.visualization import astropy_mpl_style
    import astropy.units as u
    import numpy as np
    import matplotlib.pyplot as plt
    plt.style.use(astropy_mpl_style)
    import gala.potential as gp
    from gala.units import dimensionless

    pot = gp.PlummerPotential(m=1., b=1., units=dimensionless)

    def sample_r(size=1):
        mu = np.random.random(size=size)
        return 1 / np.sqrt(mu**(-2/3) - 1)

    n_samples = 16384
    r = sample_r(size=n_samples)

    bins = np.logspace(-2, 3, 128)
    bin_cen = (bins[1:] + bins[:-1]) / 2.
    H,edges = np.histogram(r, bins=bins, weights=np.zeros_like(r) + pot.parameters['m']/r.size)

    V = 4/3.*np.pi*(bins[1:]**3 - bins[:-1]**3)

    _r = np.logspace(-2, 2, 1024)
    q = np.zeros((3,_r.size))
    q[0] = _r

    fig = plt.figure(figsize=(6,4))
    plt.loglog(_r, pot.density(q), marker=None, label='True profile', color='#cccccc', lw=3)
    plt.loglog(bin_cen, H / V, marker=None, label='Particles', color='k')
    plt.legend(loc='lower left')
    plt.xlim(1E-2, 1E2)
    plt.xlabel('$r$')
    plt.ylabel(r'$\rho(r)$')
    fig.tight_layout()

With the above, we now have sampled spherical radii that follow the desired
density profile. To compute the expansion coefficients needed to represent this
density using SCF with Hernquist radial functions, we first need to convert to
3D cartesian positions. We will distribute these particles uniformly in angles::

    φ = np.random.uniform(0, 2*np.pi, size=n_samples)
    θ = np.arccos(2*np.random.random(size=n_samples) - 1)

    xyz = np.zeros((n_samples, 3))
    xyz[:,0] = r * np.cos(φ) * np.sin(θ)
    xyz[:,1] = r * np.sin(φ) * np.sin(θ)
    xyz[:,2] = r * np.cos(θ)

    plt.figure(figsize=(5,5))
    plt.plot(xyz[:,0], xyz[:,1], linestyle='none',
             marker=',', alpha=0.25, color='k')
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.xlabel('$x$')
    plt.ylabel('$y$')

.. plot::
    :align: center
    :context: close-figs

    φ = np.random.uniform(0, 2*np.pi, size=n_samples)
    θ = np.arccos(2*np.random.random(size=n_samples) - 1)

    xyz = np.zeros((n_samples, 3))
    xyz[:,0] = r * np.cos(φ) * np.sin(θ)
    xyz[:,1] = r * np.sin(φ) * np.sin(θ)
    xyz[:,2] = r * np.cos(θ)

    plt.figure(figsize=(5,5))
    plt.plot(xyz[:,0], xyz[:,1], linestyle='none',
             marker=',', alpha=0.25, color='k')
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.xlabel('$x$')
    plt.ylabel('$y$')

To compute the expansion coefficients, we then pass the positions ``xyz`` and
masses of each "particle" to `~biff.scf.compute_coeffs_discrete`. We will
generate an array of masses that sum to 1, per our choice of units above. To
start, we'll assume that the particle distribution has spherical symmetry and
ignore terms with :math:`l>0`. We'll then plot the magnitude of the coefficients
as a function of :math:`n` (but we'll ignore the sine terms, :math:`T_{nlm}` for
this example)::

    mass = np.ones(n_samples) / n_samples
    S,T = bscf.compute_coeffs_discrete(xyz, mass=mass, nmax=16, lmax=0, r_s=1.)

    plt.semilogy(np.abs(S[:,0,0]), marker=None, lw=2)
    plt.xlabel("$n$")
    plt.ylabel("$S_{n00}$")
    plt.tight_layout()

.. plot::
    :align: center
    :context: close-figs

    import biff.scf as bscf

    mass = np.ones(n_samples) / n_samples
    S,T = bscf.compute_coeffs_discrete(xyz, mass=mass, nmax=20, lmax=0, r_s=1.)

    plt.figure(figsize=(6,4))
    plt.semilogy(np.abs(S[:,0,0]), marker=None, lw=2)
    plt.xlabel("$n$")
    plt.ylabel("$S_{n00}$")
    plt.tight_layout()

In addition to computing the coefficient values, we can also compute the
variances of the coefficients. This will let us estimate the signal-to-noise of
each expansion term and will aid us in deciding when to truncate the expansion
(see [W96]_ for the methodology and reasoning behind this)::

    (S,varS),(T,varT) = bscf.compute_coeffs_discrete(xyz, mass=mass, r_s=1.,
                                                     nmax=20, lmax=0,
                                                     compute_var=True)

    signal_to_noise = np.sqrt(S**2 / varS)

    plt.plot(signal_to_noise, marker=None, lw=2)
    plt.xlabel("$n$")
    plt.ylabel("$S/N$")

.. plot::
    :align: center
    :context: close-figs

    (S,varS),(T,varT) = bscf.compute_coeffs_discrete(xyz, mass=mass, r_s=1.,
                                                     nmax=20, lmax=0,
                                                     compute_var=True)

    signal_to_noise = np.sqrt(S[:,0,0]**2 / varS[:,0,0])

    plt.figure(figsize=(6,4))
    plt.plot(signal_to_noise, marker=None, lw=2)
    plt.xlabel("$n$")
    plt.ylabel("$S/N$")
    plt.tight_layout()

.. _coeff-analytic:

Computing expansion coefficients for an analytic density
--------------------------------------------------------

Here's another example of integrating the
`Lorenz equations <https://en.wikipedia.org/wiki/Lorenz_system>`_, a 3D
nonlinear system::

    >>> print("TODO")
    TODO

.. plot::
    :align: center

    import astropy.units as u
    import matplotlib.pyplot as pl
    import numpy as np
    import gala.integrate as gi

References
----------
.. [W96] http://dx.doi.org/10.1086/177902
.. [HMV11] http://www.artcompsci.org/kali/vol/plummer/volume11.pdf
