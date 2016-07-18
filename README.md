Biff
====

Basis Function Expansions of gravitational potentials.

[![Coverage Status](https://coveralls.io/repos/github/adrn/biff/badge.svg?branch=master)](https://coveralls.io/github/adrn/biff?branch=master)
[![Build status](http://img.shields.io/travis/adrn/biff/master.svg?style=flat)](http://travis-ci.org/adrn/biff)
[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/adrn/biff/blob/master/LICENSE)
[![Documentation Status](https://readthedocs.org/projects/gala-astro/badge/?version=latest)](http://gala-astro.readthedocs.io/en/latest/?badge=latest)
<!-- [![PyPI](https://badge.fury.io/py/biff.svg)](https://badge.fury.io/py/biff) -->

Documentation
=============

<!-- [Read the docs](http://gala.adrian.pw) -->
Working on it! Check back soon

Why?
====

Dynamical modeling of the Galactic halo is often done with simple, analytic potential forms; these are unrealistic and artificially restrictive models. If we want to understand the degeneracies and inferential power of the data at hand, we need to move to more flexible models that incorporate our ignorance on the shape of the Galactic potential at distances >~20 kpc. Representing the potential using a basis function expansion is one way to do this.

Dependencies
============

To build the Cython code, you'll need the [GNU scientific library](http://www.gnu.org/software/gsl/). On a Mac, I recommend installing this with [Homebrew](http://brew.sh/). Then, you should be able to just do:

```bash
python setup.py install
```
