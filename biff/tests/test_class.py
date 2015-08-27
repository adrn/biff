# coding: utf-8

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import astropy.units as u
from astropy.constants import G as _G
G = _G.decompose([u.kpc,u.Myr,u.Msun]).value
import numpy as np
import gary.potential as gp
from gary.units import galactic

# Project
from .._bfe import SCFPotential

def test_hernquist():
    nmax = 6
    lmax = 2

    M = 1E10
    r_s = 3.5

    cos_coeff = np.zeros((nmax+1,lmax+1,lmax+1))
    sin_coeff = np.zeros((nmax+1,lmax+1,lmax+1))
    cos_coeff[0,0,0] = 1.
    scf_potential = SCFPotential(m=M, r_s=r_s,
                                 cos_coeff=cos_coeff, sin_coeff=sin_coeff,
                                 units=galactic)

    nbins = 128
    rr = np.linspace(0.1,10.,nbins)
    xyz = np.zeros((nbins,3))
    xyz[:,0] = rr * np.cos(np.pi/4.) * np.sin(np.pi/4.)
    xyz[:,1] = rr * np.sin(np.pi/4.) * np.sin(np.pi/4.)
    xyz[:,2] = rr * np.cos(np.pi/4.)

    hernquist = gp.HernquistPotential(m=M, c=r_s, units=galactic)

    bfe_pot = scf_potential.value(xyz)
    true_pot = hernquist.value(xyz)
    np.testing.assert_allclose(bfe_pot, true_pot)

    bfe_grad = scf_potential.gradient(xyz)
    true_grad = hernquist.gradient(xyz)
    np.testing.assert_allclose(bfe_grad, true_grad)
