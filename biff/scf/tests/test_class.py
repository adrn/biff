# coding: utf-8

from __future__ import division, print_function

# Third-party
import astropy.units as u
from astropy.constants import G as _G
G = _G.decompose([u.kpc,u.Myr,u.Msun]).value
import numpy as np
import gala.potential as gp
from gala.units import galactic
from gala.potential.potential.tests.helpers import PotentialTestBase
from gala.potential.potential.io import load

# Project
from .. import _bfe

def test_hernquist():
    nmax = 6
    lmax = 2

    M = 1E10
    r_s = 3.5

    cos_coeff = np.zeros((nmax+1,lmax+1,lmax+1))
    sin_coeff = np.zeros((nmax+1,lmax+1,lmax+1))
    cos_coeff[0,0,0] = 1.
    scf_potential = _bfe.SCFPotential(m=M, r_s=r_s,
                                      Snlm=cos_coeff, Tnlm=sin_coeff,
                                      units=galactic)

    q = np.array([[1.,2.,3.]])
    t = np.array([0.])
    print("-----------------")
    print(scf_potential.c_instance.energy(q, t))

    return

    # nbins = 128
    # rr = np.linspace(0.1,10.,nbins)
    # xyz = np.zeros((3,nbins))
    # xyz[0] = rr * np.cos(np.pi/4.) * np.sin(np.pi/4.)
    # xyz[1] = rr * np.sin(np.pi/4.) * np.sin(np.pi/4.)
    # xyz[2] = rr * np.cos(np.pi/4.)

    # HACK
    xyz = (np.array([1.,2.,3.]) * u.kpc).reshape(3,1)
    scf_potential._energy(np.ascontiguousarray(xyz.value.T), t=np.array([0.]))
    # scf_potential.gradient(xyz)
    return

    hernquist = gp.HernquistPotential(m=M, c=r_s, units=galactic)

    bfe_pot = scf_potential.energy(xyz).value
    true_pot = hernquist.energy(xyz).value
    np.testing.assert_allclose(bfe_pot, true_pot)

    # print(scf_potential.gradient(xyz).value)
    return

    bfe_grad = scf_potential.gradient(xyz).value
    true_grad = hernquist.gradient(xyz).value
    np.testing.assert_allclose(bfe_grad, true_grad)

# class TestSCFPotential(PotentialTestBase):
#     nmax = 6
#     lmax = 2
#     Snlm = np.zeros((nmax+1,lmax+1,lmax+1))
#     Tnlm = np.zeros((nmax+1,lmax+1,lmax+1))
#     Snlm[0,0,0] = 1.
#     Snlm[2,0,0] = 0.5
#     Snlm[4,0,0] = 0.25

#     potential = _bfe.SCFPotential(m=1E11*u.Msun, r_s=10*u.kpc,
#                                   Snlm=Snlm, Tnlm=Tnlm, units=galactic)
#     w0 = [4.0,0.7,-0.9,0.0352238,0.1579493,0.02]

#     def test_save_load(self, tmpdir):
#         fn = str(tmpdir.join("{}.yml".format(self.name)))
#         self.potential.save(fn)
#         p = load(fn, module=_bfe)
#         p.value(self.w0[:self.w0.size//2])
