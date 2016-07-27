Biff
====

Basis function expansions of gravitational potentials.

[![Coverage Status](https://coveralls.io/repos/github/adrn/biff/badge.svg?branch=master)](https://coveralls.io/github/adrn/biff?branch=master)
[![Build status](http://img.shields.io/travis/adrn/biff/master.svg?style=flat)](http://travis-ci.org/adrn/biff)
[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/adrn/biff/blob/master/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/gala-astro/badge/?version=latest)](http://gala-astro.readthedocs.io/en/latest/?badge=latest)
[![PyPI](https://badge.fury.io/py/biff.svg)](https://badge.fury.io/py/biff)

Documentation
=============

[Read the docs.](http://biff.readthedocs.io)

Installation dependencies
=========================

The internals of Biff are implemented in C and Cython. To build the Cython code, you'll need to
have installed [GNU scientific library (GSL)](http://www.gnu.org/software/gsl/). On a Mac, I
recommend installing this with [Homebrew](http://brew.sh/) or [anaconda](http://anaconda.org). With
anaconda, you can install GSL with:

    conda -c https://anaconda.org/asmeurer/gsl install gsl

You'll also need the following packages:

    - [Cython](https://github.com/cython/cython)
    - [Numpy](https://github.com/numpy/numpy)
    - [Scipy](https://github.com/scipy/scipy)
    - [Astropy](https://github.com/astropy/astropy)
    - [Gala](https://github.com/adrn/gala) (install with `pip install astro-gala`)

Installation
============

The easiest way to install the latest stable version of Biff is via pip:

    pip install biff

To get the latest development version, you can clone this repository and do:

    python setup.py install

License
=======

[See the license.](https://github.com/adrn/biff/blob/master/LICENSE)
