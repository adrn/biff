# Licensed under a 3-clause BSD style license - see PYFITS.rst
from __future__ import absolute_import

import sys
from os import path
from distutils.core import Extension
from distutils.sysconfig import get_python_inc
from astropy_helpers import setup_helpers

def get_extensions():
    exts = []

    include_dirs = []

    # C header files in this project
    include_dirs.append('biff/scf/src')

    # malloc
    mac_incl_path = "/usr/include/malloc"
    include_dirs.append(mac_incl_path)

    # we use the string 'numpy' because astropy_helpers will auto-fill
    include_dirs.append('numpy')

    try:
        # Get gala path -- for egg_info build
        import gala
        gala_base_path = path.split(gala.__file__)[0]
        gala_path = path.join(gala_base_path, 'potential')
    except ImportError:
        gala_path = None

    if gala_path:
        include_dirs.append(gala_path)

    # Some READTHEDOCS hacks - see
    # https://github.com/pyFFTW/pyFFTW/pull/161/files
    # https://github.com/pyFFTW/pyFFTW/pull/162/files
    include_dirs.append(path.join(sys.prefix, 'include'))
    library_dirs = path.join(sys.prefix, 'lib')

    # Python include / lib
    include_dirs.append(path.abspath(path.join(get_python_inc(),"..")))

    coeff_cfg = setup_helpers.DistutilsExtensionArgs()
    coeff_cfg['include_dirs'] = coeff_cfg['include_dirs'] + include_dirs
    coeff_cfg['library_dirs'] = coeff_cfg['library_dirs'] + library_dirs
    coeff_cfg['sources'].append('biff/scf/computecoeff.pyx')
    coeff_cfg['sources'].append('biff/scf/src/bfe_helper.c')
    # coeff_cfg['sources'].append('biff/src/coeff_helper.c')
    coeff_cfg['libraries'] = ['gsl', 'gslcblas']
    coeff_cfg['extra_compile_args'] = ['--std=gnu99']
    exts.append(Extension('biff.scf._computecoeff', **coeff_cfg))

    bfe_cfg = setup_helpers.DistutilsExtensionArgs()
    bfe_cfg['include_dirs'] = bfe_cfg['include_dirs'] + include_dirs
    bfe_cfg['library_dirs'] = bfe_cfg['library_dirs'] + library_dirs

    if gala_path is not None:
        bfe_cfg['include_dirs'].append(gala_path)
        bfe_cfg['sources'].append(path.join(gala_path, 'potential',
                                            'src', 'cpotential.c'))
        bfe_cfg['sources'].append(path.join(gala_path, 'potential',
                                            'builtin', 'builtin_potentials.c'))

    bfe_cfg['sources'].append('biff/scf/bfe.pyx')
    bfe_cfg['sources'].append('biff/scf/src/bfe.c')
    bfe_cfg['sources'].append('biff/scf/src/bfe_helper.c')
    bfe_cfg['libraries'] = ['gsl', 'gslcblas']
    bfe_cfg['extra_compile_args'] = ['--std=gnu99']
    exts.append(Extension('biff.scf._bfe', **bfe_cfg))

    # SCFPotential class
    cfg = setup_helpers.DistutilsExtensionArgs()
    cfg['include_dirs'] = cfg['include_dirs'] + include_dirs
    cfg['library_dirs'] = cfg['library_dirs'] + library_dirs
    cfg['extra_compile_args'].append('--std=gnu99')
    cfg['sources'].append('biff/scf/bfe_class.pyx')
    cfg['sources'].append('biff/scf/src/bfe.c')
    cfg['sources'].append('biff/scf/src/bfe_helper.c')
    cfg['libraries'] = ['gsl', 'gslcblas']
    exts.append(Extension('biff.scf.bfe_class', **cfg))

    return exts

def get_package_data():
    return {'biff.scf': ['*.pyx',
                         'tests/data/*',
                         'tests/data/*.csv',
                         'tests/data/*.dat.gz',
                         'tests/data/*.coeff',
                         '*.h', '*.pyx', '*.pxd',
                         'src/*.c', 'src/*.h']}
