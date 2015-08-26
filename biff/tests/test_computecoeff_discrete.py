# coding: utf-8

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

import os

# Third-party
import numpy as np

# Project
from ..core import compute_coeffs_discrete
from .._bfe import potential

# HACK!
test_data_path = "/Users/adrian/projects/biff/test-data/plummer.txt"

def test_plummer():
    scfbi = scfbi = np.loadtxt(test_data_path)
    m_k = scfbi[:,0]
    xyz = scfbi[:,1:4]

    G = 1.
    r_s = 1.
    M = 1.

    nmax = 10
    lmax = 0

    nlms = []
    for n in range(nmax+1):
        for l in range(lmax+1):
            for m in range(l+1):
                nlms.append([n,l,m])

    Snlm, Tnlm = np.zeros((2,nmax+1,lmax+1,lmax+1))

    for nlm in nlms:
        n,l,m = nlm
        Snlm[n,l,m],Tnlm[n,l,m] = compute_coeffs_discrete(xyz, m_k, nlm, r_s=r_s)

    # for first 100 particles
    for ix in range(100):
        ixs = np.arange(m_k.size).astype(int)
        ixs = np.delete(ixs, ix)

        val = -np.sum(m_k[ixs] / np.linalg.norm(xyz[ixs] - xyz[ix], axis=-1))
        biff_val = potential(np.ascontiguousarray(xyz[ix][None]), G, M, r_s, Snlm, Tnlm, nmax, lmax)[0]

        np.testing.assert_allclose(val, biff_val, rtol=0.05) # 5% tolerance