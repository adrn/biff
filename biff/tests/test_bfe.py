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
    xyz = np.atleast_2d(xyz)
    r = np.sqrt(np.sum(xyz**2, axis=-1))
    return M/(2*np.pi) * r_s / (r * (r+r_s)**3)

def hernquist_potential(xyz, M, r_s):
    xyz = np.atleast_2d(xyz)
    r = np.sqrt(np.sum(xyz**2, axis=-1))
    return -G*M / (r + r_s)

def hernquist_gradient(xyz, M, r_s):
    import gary.potential as gp
    p = gp.HernquistPotential(m=M, c=r_s, units=[u.kpc,u.Myr,u.Msun,u.radian])
    return p.gradient(xyz)

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
    true_dens = hernquist_density(xyz, M, r_s)
    np.testing.assert_allclose(bfe_dens, true_dens)

    bfe_pot = potential(xyz, Snlm, Tnlm, nmax, lmax, G=G, M=M, r_s=r_s)
    true_pot = hernquist_potential(xyz, M, r_s)
    np.testing.assert_allclose(bfe_pot, true_pot)

    bfe_grad = gradient(xyz, Snlm, Tnlm, nmax, lmax, G=G, M=M, r_s=r_s)
    true_grad = hernquist_gradient(xyz, M, r_s)
    np.testing.assert_allclose(bfe_grad, true_grad)
