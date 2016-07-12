# Licensed under a 3-clause BSD style license - see PYFITS.rst
from __future__ import absolute_import

import os
from distutils.core import Extension
from astropy_helpers import setup_helpers

def get_extensions():

    try:
        # Get gala path
        import gala
        gala_base_path = os.path.split(gala.__file__)[0]
        gala_path = os.path.join(gala_base_path, 'potential')
    except ImportError:
        gala_path = None

    coeff_cfg = setup_helpers.DistutilsExtensionArgs()
    coeff_cfg['include_dirs'].append('numpy')
    coeff_cfg['include_dirs'].append('biff/src')
    coeff_cfg['sources'].append('biff/computecoeff.pyx')
    coeff_cfg['sources'].append('biff/src/bfe_helper.c')
    # coeff_cfg['sources'].append('biff/src/coeff_helper.c')
    coeff_cfg['libraries'] = ['gsl', 'gslcblas']
    coeff_cfg['extra_compile_args'] = ['--std=gnu99']

    bfe_cfg = setup_helpers.DistutilsExtensionArgs()
    bfe_cfg['include_dirs'].append('numpy')
    if gala_path is not None:
        bfe_cfg['include_dirs'].append(gala_path)
    bfe_cfg['include_dirs'].append('biff/src')
    bfe_cfg['sources'].append('biff/bfe.pyx')
    bfe_cfg['sources'].append('biff/src/bfe.c')
    bfe_cfg['sources'].append('biff/src/bfe_helper.c')
    # bfe_cfg['sources'].append(os.path.join(gala_path, 'src', 'cpotential.c'))
    bfe_cfg['libraries'] = ['gsl', 'gslcblas']
    bfe_cfg['extra_compile_args'] = ['--std=gnu99']

    return [Extension('biff._computecoeff', **coeff_cfg),
            Extension('biff._bfe', **bfe_cfg)]

def get_package_data():
    return {'biff': ['*.pyx', '*/*.pyx',
                     'src/*.h', 'data/*.dat.gz', 'data/*.coeff']}
