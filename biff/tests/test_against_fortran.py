# coding: utf-8

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import os

# Third-party
import astropy.units as u
from astropy.constants import G as _G
from astropy.utils.data import get_pkg_data_filename
from astropy.tests.helper import pytest
G = _G.decompose([u.kpc,u.Myr,u.Msun]).value
import numpy as np

# Project
from .._bfe import density, potential, gradient

@pytest.mark.parametrize("basename", [
    'simple-hernquist', 'random', 'wang-zhao',
])
@pytest.mark.skipif('True')
def test_density(basename):
    pos_path = os.path.abspath(get_pkg_data_filename('../data/positions.dat.gz'))
    coeff_path = os.path.abspath(get_pkg_data_filename('../data/{0}.coeff'.format(basename)))
    accp_path = os.path.abspath(get_pkg_data_filename('../data/{0}-accp.dat.gz'.format(basename)))

    xyz = np.loadtxt(pos_path, skiprows=1)
    coeff = np.atleast_2d(np.loadtxt(coeff_path, skiprows=1))

    nmax = coeff[:,0].astype(int).max()
    lmax = coeff[:,1].astype(int).max()

    cos_coeff = np.zeros((nmax+1,lmax+1,lmax+1))
    sin_coeff = np.zeros((nmax+1,lmax+1,lmax+1))
    for row in coeff:
        n,l,m,cc,sc = row
        cos_coeff[int(n),int(l),int(m)] = cc
        sin_coeff[int(n),int(l),int(m)] = sc

    dens = density(xyz, M=1., r_s=1.,
                   Snlm=cos_coeff, Tnlm=sin_coeff,
                   nmax=nmax, lmax=lmax)

    # TODO: nothing to compare this to....

@pytest.mark.parametrize("basename", [
    'simple-hernquist', 'random', 'wang-zhao',
])
def test_potential(basename):
    pos_path = os.path.abspath(get_pkg_data_filename('../data/positions.dat.gz'))
    coeff_path = os.path.abspath(get_pkg_data_filename('../data/{0}.coeff'.format(basename)))
    accp_path = os.path.abspath(get_pkg_data_filename('../data/{0}-accp.dat.gz'.format(basename)))

    xyz = np.loadtxt(pos_path, skiprows=1)
    coeff = np.atleast_2d(np.loadtxt(coeff_path, skiprows=1))
    accp = np.loadtxt(accp_path)

    nmax = coeff[:,0].astype(int).max()
    lmax = coeff[:,1].astype(int).max()

    cos_coeff = np.zeros((nmax+1,lmax+1,lmax+1))
    sin_coeff = np.zeros((nmax+1,lmax+1,lmax+1))
    for row in coeff:
        n,l,m,cc,sc = row
        cos_coeff[int(n),int(l),int(m)] = cc
        sin_coeff[int(n),int(l),int(m)] = sc

    potv = potential(xyz, G=1., M=1., r_s=1.,
                     Snlm=cos_coeff, Tnlm=sin_coeff,
                     nmax=nmax, lmax=lmax)

    # for some reason, SCF potential is -potential
    scf_potv = -accp[:,-1]
    np.testing.assert_allclose(potv, scf_potv, rtol=1E-6)

@pytest.mark.parametrize("basename", [
    'simple-hernquist', 'random', 'wang-zhao',
])
def test_gradient(basename):
    pos_path = os.path.abspath(get_pkg_data_filename('../data/positions.dat.gz'))
    coeff_path = os.path.abspath(get_pkg_data_filename('../data/{0}.coeff'.format(basename)))
    accp_path = os.path.abspath(get_pkg_data_filename('../data/{0}-accp.dat.gz'.format(basename)))

    xyz = np.loadtxt(pos_path, skiprows=1)
    coeff = np.atleast_2d(np.loadtxt(coeff_path, skiprows=1))
    accp = np.loadtxt(accp_path)

    nmax = coeff[:,0].astype(int).max()
    lmax = coeff[:,1].astype(int).max()

    cos_coeff = np.zeros((nmax+1,lmax+1,lmax+1))
    sin_coeff = np.zeros((nmax+1,lmax+1,lmax+1))
    for row in coeff:
        n,l,m,cc,sc = row
        cos_coeff[int(n),int(l),int(m)] = cc
        sin_coeff[int(n),int(l),int(m)] = sc

    grad = gradient(xyz, G=1., M=1., r_s=1.,
                    Snlm=cos_coeff, Tnlm=sin_coeff,
                    nmax=nmax, lmax=lmax)

    # I output the acceleration from SCF when I make the files
    #   so I have no idea why I don't need a minus sign here...
    scf_grad = accp[:,:3]
    np.testing.assert_allclose(grad, scf_grad, rtol=1E-6)
