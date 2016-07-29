============
Installation
============

The easiest way to install the latest stable version of Biff is via pip::

    pip install biff

If this doesn't work, you can try to build the latest developer version as
described below.

Cloning, Building, Installing
=============================

The latest development version of Biff can be cloned from
`GitHub <https://github.com/>`_ using ``git``::

   git clone git://github.com/adrn/biff.git

To build the project (from the root of the source tree, e.g., inside
the cloned ``biff`` directory)::

    python setup.py build

To install the project::

    python setup.py install

Dependencies
============

The internals of Biff are implemented in C and Cython. To build the Cython
code, you'll need to have installed `GNU scientific library (GSL)
<http://www.gnu.org/software/gsl/>`_. On a Mac, I recommend installing this
with `Homebrew <http://brew.sh/>`_ or `anaconda <http://anaconda.org>`_. With
anaconda, you can install GSL with::

    conda -c https://anaconda.org/asmeurer/gsl install gsl

You'll also need the following packages:

- `Cython <https://github.com/cython/cython>`_
- `Numpy <https://github.com/numpy/numpy>`_
- `Scipy <https://github.com/scipy/scipy>`_
- `Astropy`_
- `Gala <https://github.com/adrn/gala>`_ (install with ``pip install astro-gala``)
