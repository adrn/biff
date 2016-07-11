# coding: utf-8

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import astropy.units as u
from astropy.constants import G as _G
G = _G.decompose([u.kpc,u.Myr,u.Msun]).value
import numpy as np

# Project
from .._bfe import density, potential, gradient

# Check that we get A000=1. for putting in hernquist density
def hernquist_density(xyz, M, r_s):
    r = np.sqrt(np.sum(xyz**2, axis=0))
    return M/(2*np.pi) * r_s / (r * (r+r_s)**3)

def hernquist_potential(xyz, M, r_s):
    r = np.sqrt(np.sum(xyz**2, axis=0))
    return -G*M / (r + r_s)

def hernquist_gradient(xyz, M, r_s):
    import gala.potential as gp
    p = gp.HernquistPotential(m=M, c=r_s, units=[u.kpc,u.Myr,u.Msun,u.radian])
    return p.gradient(xyz).value

def test_hernquist():
    nmax = 6
    lmax = 2

    Snlm = np.zeros((nmax+1,lmax+1,lmax+1))
    Tnlm = np.zeros((nmax+1,lmax+1,lmax+1))
    Snlm[0,0,0] = 1.

    M = 1E10
    r_s = 3.5

    nbins = 128
    rr = np.linspace(0.1,10.,nbins)
    xyz = np.zeros((nbins,3))
    xyz[:,0] = rr * np.cos(np.pi/4.) * np.sin(np.pi/4.)
    xyz[:,1] = rr * np.sin(np.pi/4.) * np.sin(np.pi/4.)
    xyz[:,2] = rr * np.cos(np.pi/4.)

    bfe_dens = density(xyz, Snlm, Tnlm, nmax, lmax, M=M, r_s=r_s)
    true_dens = hernquist_density(xyz.T, M, r_s)
    np.testing.assert_allclose(bfe_dens, true_dens)

    bfe_pot = potential(xyz, Snlm, Tnlm, nmax, lmax, G=G, M=M, r_s=r_s)
    true_pot = hernquist_potential(xyz.T, M, r_s)
    np.testing.assert_allclose(bfe_pot, true_pot)

    bfe_grad = gradient(xyz, Snlm, Tnlm, nmax, lmax, G=G, M=M, r_s=r_s)
    true_grad = hernquist_gradient(xyz.T, M, r_s)
    np.testing.assert_allclose(bfe_grad.T, true_grad)

# ----------------------------------------------------------------------------

def pure_py(xyz, Snlm, Tnlm, nmax, lmax):
    from scipy.special import lpmv, gegenbauer
    from math import factorial as f

    Plm = lambda l,m,costh: lpmv(m, l, costh)
    Ylmth = lambda l,m,costh: np.sqrt((2*l+1)/(4 * np.pi) * f(l-m)/f(l+m)) * Plm(l,m,costh)

    twopi = 2*np.pi
    sqrt4pi = np.sqrt(4*np.pi)

    r = np.sqrt(np.sum(xyz**2, axis=0))
    X = xyz[2]/r # cos(theta)
    phi = np.arctan2(xyz[1], xyz[0])
    xsi = (r - 1) / (r + 1)

    density = 0
    potenti = 0
    for l in range(lmax+1):
        r_term1 = r**l / (r*(1+r)**(2*l+3))
        r_term2 = r**l / (1+r)**(2*l+1)
        for m in range(l+1):
            for n in range(nmax+1):
                Cn = gegenbauer(n, 2*l+3/2)
                Knl = 0.5 * n * (n+4*l+3) + (l+1)*(2*l+1)
                rho_nl = Knl / twopi * sqrt4pi * r_term1 * Cn(xsi)
                phi_nl = -sqrt4pi * r_term2 * Cn(xsi)

                density += rho_nl * Ylmth(l,m,X) * (Snlm[n,l,m]*np.cos(m*phi) +
                                                    Tnlm[n,l,m]*np.sin(m*phi))
                potenti += phi_nl * Ylmth(l,m,X) * (Snlm[n,l,m]*np.cos(m*phi) +
                                                    Tnlm[n,l,m]*np.sin(m*phi))

    return density, potenti

def test_pure_py():

    nmax = 6
    lmax = 4
    xyz = np.ascontiguousarray(np.array([[1.,0.,1.], [1.,1.,0.], [0.,1.,1.]]).T)

    # Snlm = np.random.uniform(-1,1,size=(nmax+1,lmax+1,lmax+1))
    Snlm = np.zeros((nmax+1,lmax+1,lmax+1))
    Snlm[0,0,0] = 1.
    Snlm[2,0,0] = 1E-3
    Snlm[4,0,0] = 1E-5
    Snlm[6,0,0] = 1E-6
    Tnlm = np.zeros_like(Snlm)

    py_den,py_pot = pure_py(xyz, Snlm, Tnlm, nmax, lmax)

    cy_den = density(xyz, Snlm, Tnlm, nmax, lmax, M=1., r_s=1.)
    cy_pot = potential(xyz, Snlm, Tnlm, nmax, lmax, G=1., M=1., r_s=1.)

    print("Density:")
    print("\tPython", py_den)
    print("\tCython", cy_den)
    print("-"*64)
    print("Potential:")
    print("\tPython", py_pot)
    print("\tCython", cy_pot)
