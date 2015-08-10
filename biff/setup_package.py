# Licensed under a 3-clause BSD style license - see PYFITS.rst
from __future__ import absolute_import

from distutils.core import Extension

from astropy_helpers import setup_helpers


def get_extensions():
    cfg = setup_helpers.DistutilsExtensionArgs()
    # 'numpy' will be replaced with the proper path to the numpy includes
    cfg['include_dirs'].append('numpy')
    cfg['include_dirs'].append('biff/src')
    cfg['sources'].append('biff/computecoeff.pyx')
    cfg['libraries'] = ['gsl', 'gslcblas']

    return [Extension('biff._bfe', **cfg)]
