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
    # cfg['sources'].append('cextern/brent.c')
    # cfg['sources'].append('cextern/simpson.c')
    ext = Extension('biff._bfe', **cfg)

    return [ext]

# are the c files provided by some external library?  If so, should add it here
# so that users have the option of using that instead of the builtin one
def get_external_libraries():
    return ['gsl', 'gslcblas']
