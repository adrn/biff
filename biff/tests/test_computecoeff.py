# coding: utf-8

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import astropy.units as u
from astropy.constants import G as _G
G = _G.decompose([u.kpc,u.Myr,u.Msun]).value
import numpy as np
import gala.potential as gp
from gala.units import galactic

# Project
from ..core import compute_coeffs
from .._bfe import potential, density, gradient

# Check that we get A000=1. for putting in hernquist density
def hernquist_density(x, y, z, M, r_s):
    r = np.sqrt(x**2 + y**2 + z**2)
    return M/(2*np.pi) * r_s / (r * (r+r_s)**3)

def test_hernquist():
    for M in [1E5, 1E10]:
        for r_s in np.logspace(-1,2,4):
            (S,Serr),(T,Terr) = compute_coeffs(hernquist_density, nlm=[0,0,0], M=M, r_s=r_s, args=(M, r_s))

            np.testing.assert_allclose(S, 1.)
            np.testing.assert_allclose(Serr, 0., atol=1E-10)

            np.testing.assert_allclose(T, 0.)
            np.testing.assert_allclose(Terr, 0., atol=1E-10)

# ----------------------------------------------------------------------------

def _plummer_density(x, y, z, M, r_s):
    r2 = x*x + y*y + z*z
    return (3*M / (4*np.pi*r_s**3)) * (1 + r2/r_s**2)**(-5/2.)

def test_plummer():
    true_M = 1/G
    true_r_s = 1.

    x = np.logspace(-2,1,512)
    xyz = np.zeros((len(x),3))
    xyz[:,0] = x

    pot = gp.PlummerPotential(m=true_M, b=true_r_s, units=galactic)
    true_pot = pot.value(xyz.T).value
    true_dens = pot.density(xyz.T).value
    true_grad = pot.gradient(xyz.T).value.T

    nmax = 16
    lmax = 0

    Snlm = np.zeros((nmax+1,lmax+1,lmax+1))
    Serr = np.zeros((nmax+1,lmax+1,lmax+1))
    Tnlm = np.zeros((nmax+1,lmax+1,lmax+1))
    Terr = np.zeros((nmax+1,lmax+1,lmax+1))

    nlms = []
    for n in range(nmax+1):
        for l in range(lmax+1):
            for m in range(l+1):
                nlms.append([n,l,m])

    for nlm in nlms:
        n,l,m = nlm
        print(n,l,m)
        (S,S_err),(T,T_err) = compute_coeffs(_plummer_density, nlm=nlm,
                                             M=true_M, r_s=true_r_s, args=(true_M,true_r_s),
                                             epsrel=1E-9)
        Snlm[n,l,m] = S
        Serr[n,l,m] = S_err
        Tnlm[n,l,m] = T
        Terr[n,l,m] = T_err

    bfe_dens = density(xyz, Snlm, Tnlm, nmax, lmax, true_M, true_r_s)
    bfe_pot = potential(xyz, Snlm, Tnlm, nmax, lmax, G, true_M, true_r_s)
    bfe_grad = gradient(xyz, Snlm, Tnlm, nmax, lmax, G, true_M, true_r_s)

    assert np.allclose(true_dens, bfe_dens, rtol=2E-3)
    assert np.allclose(true_pot, bfe_pot, rtol=1E-6)
    assert np.allclose(true_grad[:,0], bfe_grad[:,0], rtol=5E-3)
    # print(np.abs((bfe_dens - true_dens) / true_dens).max())
    # print(np.abs((bfe_pot - true_pot) / true_pot).max())
    # print(np.abs((bfe_grad[:,0] - true_grad[:,0]) / true_grad[:,0]).max())

# ----------------------------------------------------------------------------

def miyamoto_nagai_density(x, y, z, M, a, b):
    R2 = x**2 + y**2
    sqrt_zb = np.sqrt(z**2+b**2)

    A = (b**2*M / (4*np.pi))
    numer = a*R2 + (a + 2*sqrt_zb)*(a + sqrt_zb)**2
    denom = (R2 + (a + sqrt_zb)**2)**2.5 * sqrt_zb**3

    return A*numer/denom

def test_miyamoto_nagai():
    M = 1E10
    a = 1.
    b = 1.

    x = np.logspace(-2,1,512)
    xyz = np.zeros((len(x),3))
    xyz[:,0] = x

    pot = gp.MiyamotoNagaiPotential(m=M, a=a, b=b, units=galactic)
    true_pot = pot.value(xyz.T).value
    true_dens = pot.density(xyz.T).value
    true_grad = pot.gradient(xyz.T).value

    nmax = 16
    lmax = 0

    Snlm = np.zeros((nmax+1,lmax+1,lmax+1))
    Snlm_e = np.zeros((nmax+1,lmax+1,lmax+1))
    Tnlm = np.zeros((nmax+1,lmax+1,lmax+1))
    Tnlm_e = np.zeros((nmax+1,lmax+1,lmax+1))
    for n in range(nmax+1):
        for l in range(lmax+1):
            for m in range(l+1):
                nlm = [n,l,m]
                (S,Serr),(T,Terr) = compute_coeffs(miyamoto_nagai_density,
                                                   nlm=nlm,
                                                   M=M, r_s=a,
                                                   args=(M,a,b))
                Snlm[n,l,m] = S
                Snlm_e[n,l,m] = Serr
                Tnlm[n,l,m] = T
                Tnlm_e[n,l,m] = Terr

    # Snlm[(np.abs(Snlm) < np.abs(Snlm_e))] = 0.
    # Tnlm[(np.abs(Tnlm) < np.abs(Tnlm_e))] = 0.

    bfe_dens = density(xyz, Snlm, Tnlm, nmax, lmax, M, a)
    bfe_pot = potential(xyz, Snlm, Tnlm, nmax, lmax, G, M, a)
    bfe_grad = gradient(xyz, Snlm, Tnlm, nmax, lmax, G, M, a)

    import matplotlib.pyplot as pl

    fig,axes = pl.subplots(3, 1, figsize=(6,12), sharex=True)

    axes[0].loglog(x, true_dens)
    axes[0].loglog(x, bfe_dens)

    axes[1].semilogx(x, true_pot)
    axes[1].semilogx(x, bfe_pot)

    axes[2].semilogx(x, true_grad[0])
    axes[2].semilogx(x, bfe_grad[:,0])

    pl.show()
