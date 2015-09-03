# coding: utf-8

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import astropy.units as u
from astropy.constants import G as _G
G = _G.decompose([u.kpc,u.Myr,u.Msun]).value
import numpy as np

# Project
from ..core import compute_coeffs
from .._bfe import potential, density

# Check that we get A000=1. for putting in hernquist density
def hernquist_density(x, y, z, args):
    M, r_s = args
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

# -------------------------------------------------------------------
# TODO: these are failing -- why do I suck?
def miyamoto_nagai_density(x, y, z, args):
    M, a, b = args
    R2 = x**2 + y**2
    sqrt_zb = np.sqrt(z**2+b**2)

    A = (b**2*M / (4*np.pi))
    numer = a*R2 + (a + 2*sqrt_zb)*(a + sqrt_zb)**2
    denom = (R2 + (a + sqrt_zb)**2)**2.5 * sqrt_zb**3

    return A*numer/denom

def miyamoto_nagai_potential(x, y, z, args):
    M, a, b = args
    return -G * M / np.sqrt(x**2 + y**2 + (a + np.sqrt(z**2 + b**2))**2)

def test_miyamoto_nagai():
    M = 1./G
    a = 1.
    b = 0.9

    nmax = 3
    lmax = 6

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

    Snlm[(np.abs(Snlm) < np.abs(Snlm_e))] = 0.
    Tnlm[(np.abs(Tnlm) < np.abs(Tnlm_e))] = 0.

    xyz = np.array([[0.7, 0., 0.1]])
    scf_potv = potential(xyz, Snlm, Tnlm, nmax=nmax, lmax=lmax, G=G, M=M, r_s=a)
    tru_potv = miyamoto_nagai_potential(xyz[:,0],xyz[:,1],xyz[:,2],(M,a,b))

    print(scf_potv)
    print(tru_potv)
