# Biff
Basis Function Expansions of gravitational potentials.

## Why?
Galactic dynamical modeling is often done with simple, analytic potential forms --- these are unrealistic and artificially restrictive models. If we want to understand the degeneracies and inferential power in our data, we need to move to more flexible models that encompass our **complete ignorance** on the shape of the Galactic potential at $r > 10~{\rm kpc}$. Representing the potential using a basis function expansion is one way to do this.

## How?
To build the Cython code, you'll need the [GNU scientific library](http://www.gnu.org/software/gsl/). On a Mac, I recommend installing this with [Homebrew](http://brew.sh/). Then, you should be able to just do:

```bash
python setup.py build_ext
python setup.py install
```

to install Biff.
