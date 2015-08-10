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

# API:
def hernquist_density(x, y, z, M, r_s, args):
    r = np.sqrt(x**2+y**2+z**2)
    return M/(2*np.pi) * r_s / (r * (r+r_s)**3)

def test():
    v = compute_Anlm(hernquist_density, nlm=[0,0,0], M=1E10, r_s=2.)
    print(v)
