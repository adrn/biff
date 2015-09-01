# Licensed under a 3-clause BSD style license - see PYFITS.rst
from __future__ import absolute_import

from distutils.core import Extension
from astropy_helpers import setup_helpers

def get_extensions():

    coeff_cfg = setup_helpers.DistutilsExtensionArgs()
    coeff_cfg['include_dirs'].append('numpy')
    coeff_cfg['include_dirs'].append('biff/src')
    coeff_cfg['sources'].append('biff/computecoeff.pyx')
    coeff_cfg['sources'].append('biff/src/bfe_helper.c')
    coeff_cfg['libraries'] = ['gsl', 'gslcblas']

    bfe_cfg = setup_helpers.DistutilsExtensionArgs()
    bfe_cfg['include_dirs'].append('numpy')
    bfe_cfg['include_dirs'].append('biff/src')
    bfe_cfg['sources'].append('biff/bfe.pyx')
    bfe_cfg['libraries'] = ['gsl', 'gslcblas']

    return [Extension('biff._computecoeff', **coeff_cfg),
            Extension('biff._bfe', **bfe_cfg)]

def get_package_data():
    return {'biff': ['src/*.h', 'data/*.dat.gz']}
