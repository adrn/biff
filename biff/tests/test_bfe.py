# coding: utf-8

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import astropy.units as u
from astropy.constants import G as _G
G = _G.decompose([u.kpc,u.Myr,u.Msun]).value
import numpy as np

# Project
from .._bfe import density, potential

# Check that we get A000=1. for putting in hernquist density
def hernquist_density(xyz, M, r_s):
    xyz = np.atleast_2d(xyz)
    r = np.sqrt(np.sum(xyz**2, axis=-1))
    return M/(2*np.pi) * r_s / (r * (r+r_s)**3)

def hernquist_potential(xyz, M, r_s):
    xyz = np.atleast_2d(xyz)
    r = np.sqrt(np.sum(xyz**2, axis=-1))
    return -G*M / (r + r_s)

def test_hernquist():
    nmax = 6
    lmax = 2

    Anlm = np.zeros((nmax+1,lmax+1,lmax+1))
    Anlm[0,0,0] = 1.

    M = 1E10
    r_s = 3.5

    nbins = 128
    rr = np.linspace(0.1,10.,nbins)
    xyz = np.zeros((nbins,3))
    xyz[:,0] = rr

    bfe_dens = density(xyz, M, r_s, Anlm, nmax, lmax)
    true_dens = hernquist_density(xyz, M, r_s)
    np.testing.assert_allclose(bfe_dens, true_dens)

    bfe_pot = potential(xyz, G, M, r_s, Anlm, nmax, lmax)
    true_pot = hernquist_potential(xyz, M, r_s)
    np.testing.assert_allclose(bfe_pot, true_pot)
