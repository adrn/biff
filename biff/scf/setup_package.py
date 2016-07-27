# Licensed under a 3-clause BSD style license - see PYFITS.rst
from __future__ import absolute_import

from os.path import join, split, abspath
from distutils.core import Extension
from distutils.sysconfig import get_python_inc
from astropy_helpers import setup_helpers

def get_extensions():

    try:
        # Get gala path -- for egg_info build
        import gala
        gala_base_path = split(gala.__file__)[0]
        gala_path = join(gala_base_path, 'potential')
    except ImportError:
        gala_path = None

    py_inc = abspath(join(get_python_inc(),".."))
    py_lib = abspath(join(get_python_inc(),"..","..","lib"))

    coeff_cfg = setup_helpers.DistutilsExtensionArgs()
    coeff_cfg['include_dirs'].append('numpy')
    coeff_cfg['include_dirs'].append('biff/scf/src')
    coeff_cfg['include_dirs'].append(py_inc) # for gsl
    # coeff_cfg['library_dirs'].append(py_lib) # for gsl
    coeff_cfg['sources'].append('biff/scf/computecoeff.pyx')
    coeff_cfg['sources'].append('biff/scf/src/bfe_helper.c')
    # coeff_cfg['sources'].append('biff/src/coeff_helper.c')
    coeff_cfg['libraries'] = ['gsl', 'gslcblas']
    coeff_cfg['extra_compile_args'] = ['--std=gnu99']

    bfe_cfg = setup_helpers.DistutilsExtensionArgs()
    bfe_cfg['include_dirs'].append('numpy')
    if gala_path is not None:
        bfe_cfg['include_dirs'].append(gala_path)
    bfe_cfg['include_dirs'].append('biff/scf/src')
    bfe_cfg['include_dirs'].append(py_inc) # for gsl
    # bfe_cfg['library_dirs'].append(py_lib) # for gsl
    bfe_cfg['sources'].append('biff/scf/bfe.pyx')
    bfe_cfg['sources'].append('biff/scf/src/bfe.c')
    bfe_cfg['sources'].append('biff/scf/src/bfe_helper.c')
    # bfe_cfg['sources'].append(os.path.join(gala_path, 'src', 'cpotential.c'))
    bfe_cfg['libraries'] = ['gsl', 'gslcblas']
    bfe_cfg['extra_compile_args'] = ['--std=gnu99']

    return [Extension('biff.scf._computecoeff', **coeff_cfg),
            Extension('biff.scf._bfe', **bfe_cfg)]

def get_package_data():
    return {'biff': ['*.pyx', '*/*.pyx',
                     'scf/src/*.h', 'scf/data/*.dat.gz', 'scf/data/*.coeff']}
