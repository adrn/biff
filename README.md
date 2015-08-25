# Biff
Basis Function Expansions of gravitational potentials.

## Why?
Dynamical modeling of the Galactic halo is often done with simple, analytic potential forms; these are unrealistic and artificially restrictive models. If we want to understand the degeneracies and inferential power of the data at hand, we need to move to more flexible models that incorporate our ignorance on the shape of the Galactic potential at distances >~20 kpc. Representing the potential using a basis function expansion is one way to do this.

## How?
To build the Cython code, you'll need the [GNU scientific library](http://www.gnu.org/software/gsl/). On a Mac, I recommend installing this with [Homebrew](http://brew.sh/). Then, you should be able to just do:

```bash
python setup.py install
```
