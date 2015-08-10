# coding: utf-8

"""
    Note: in order to run this script, you must be in the 'examples'
    directory, and you must have biff installed (e.g., using
    python setup.py install).

    This script shows how one can compute expansion coefficients
    in parallel using either Python's built-in multiprocessing
    library (for multi-core processors), or using OpenMPI.

    To run this script using multiprocessing, call the module with
    the '-n' argument to specify the number of cores to run with::

        python parallel_compute_Anlm.py -n 2

    To run this script using MPI, you must be on a machine with
    OpenMPI installed, along with mpi4py and mpipool
    (https://github.com/adrn/mpipool). Then, simply use the '--mpi'
    boolean flag combined with mpirun or mpiexec::

        mpiexec -n 4 python parallel_compute_Anlm.py --mpi

"""

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import sys
import multiprocessing
import math

# Third-party
from astropy import log as logger
import numpy as np

# Project
import biff

# Here we define a simple density function -- we'll use a Plummer model
def plummer_density(x, y, z, M, r_s, args):
    r = math.sqrt(x*x + y*y + z*z)
    return 3*M/(4*math.pi*r_s) * (1 + r**2 / r_s**2)**(-5/2.)

def worker(task):
    nlm, M, r_s = task
    return biff.compute_Anlm(plummer_density, nlm=nlm, M=M, r_s=r_s)

def main(pool):
    # maximum n and l values to use in the expansion
    nmax = 10
    lmax = 0

    # parameters for the Plummer density
    M = 1E10
    r_s = 3.5

    nlms = []
    for n in range(nmax+1):
        for l in range(lmax+1):
            for m in range(l+1):
                print(n,l,m)
                nlms.append([n,l,m])

    tasks = [[nlm, M, r_s] for nlm in nlms]
    results = pool.map(worker, tasks)
    pool.close()

    Anlm,Anlm_err = np.array(results).T
    print(Anlm)
    print(Anlm_err)

# ----------------------------------------------------------------------------
# Below here are utilities used above
# ----------------------------------------------------------------------------
class SerialPool(object):

    def close(self):
        return

    def map(self, function, tasks, callback=None):
        results = []
        for task in tasks:
            result = function(task)
            if callback is not None:
                callback(result)
            results.append(result)
        return results

def get_pool(mpi=False, threads=None, **kwargs):
    """
    Get a pool object to pass to emcee for parallel processing.
    If mpi is False and threads is None, pool is None.

    Parameters
    ----------
    mpi : bool
        Use MPI or not. If specified, ignores the threads kwarg.
    threads : int (optional)
        If mpi is False and threads is specified, use a Python
        multiprocessing pool with the specified number of threads.
    **kwargs
        Any other keyword arguments are passed through to the pool
        initializers.

    """

    if mpi:
        from mpipool import MPIPool

        # Initialize the MPI pool
        pool = MPIPool(**kwargs)

        # Make sure the thread we're running on is the master
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
        logger.debug("Running with MPI...")

    elif threads > 1:
        logger.debug("Running with multiprocessing on {} cores..."
                     .format(threads))
        pool = multiprocessing.Pool(threads, **kwargs)

    else:
        logger.debug("Running serial...")
        pool = SerialPool(**kwargs)

    return pool

if __name__ == "__main__":
    from argparse import ArgumentParser
    import logging

    # Define parser object
    parser = ArgumentParser(description="")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                        default=False, help="Be chatty! (default = False)")
    parser.add_argument("-q", "--quiet", action="store_true", dest="quiet",
                        default=False, help="Be quiet! (default = False)")

    parser.add_argument("-n", "--nthreads", dest="nthreads",
                        default=None, type=int,
                        help="Run with multiprocessing, this sets the number of threads.")
    parser.add_argument("--mpi", action="store_true", dest="mpi",
                        default=False, help="Run using MPI.")

    args = parser.parse_args()

    # Set logger level based on verbose flags
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    elif args.quiet:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

    main(pool=get_pool(mpi=args.mpi, threads=args.nthreads))
