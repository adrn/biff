# coding: utf-8

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import os
import sys

# Third-party
import numpy as np

# Project
from ..core import compute_Anlm

# Check that we get A000=1. for putting in hernquist density
def hernquist_density(x, y, z, M, r_s, args):
    r = np.sqrt(x**2+y**2+z**2)
    return M/(2*np.pi) * r_s / (r * (r+r_s)**3)

def test_hernquist():
    for M in [1E5, 1E10]:
        for r_s in np.logspace(-1,2,4):
            v,err = compute_Anlm(hernquist_density, nlm=[0,0,0], M=M, r_s=r_s)
            np.testing.assert_allclose(v, 1.)
            np.testing.assert_allclose(err, 0., atol=1E-10)
