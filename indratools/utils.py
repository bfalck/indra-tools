"""
Utility functions for the Indra suite of simulations hosted on the SciServer.

Written by Bridget Falck, 2018-2020



Inputs: 
- ``snapinput`` is an optional subset of snapshots for which redshifts and scale factors
    are desired for particle or FFT snapshots.


Methods
-------

part_snapinfo(snapinput=None)
    Reads and returns all snapnums, redshifts, and scale factors for all 64 particle 
    snapshots, or a subset of redshifts and scale factors if snapinput is set.
FFT_snapinfo(snapinput=None)
    Reads and returns all snapnums, redshifts, and scale factors for all 505 FFT 
    snapshots, or a subset of redshifts and scale factors if snapinput is set.
get_pklin()
    Reads and returns the k and P(k) values of the CAMB linear power spectrum, 
    normalized to z=0, for the Indra cosmology.

"""

import numpy as np
import pkg_resources


def part_snapinfo(snapinput=None):
    """Reads and returns all snapnums, redshifts, and scale factors for all 64 particle
    snapshots, or a subset of redshifts and scale factors if snapinput is set.
    
    Parameters
    ----------
    snapinput : int or array, optional
        Subset of snapnums for which redshifts and scale factors are desired
    
    Returns
    -------
    (snaps,) redshifts, times : ndarrays
        Arrays of redshifts and scale factors corresponding to snapnums
    """
    
    snapfile = pkg_resources.resource_filename('indratools', 'data/snapinfo.txt')
    snaps,redshifts,times = np.loadtxt(snapfile,unpack=True)
    
    if snapinput == None:
        return snaps, redshifts, times
    else:
        return redshifts[snapinput],times[snapinput]
    
    
def FFT_snapinfo(snapinput=None):
    """Reads and returns all snapnums, redshifts, and scale factors for all 505 FFT
    snapshots, or a subset of redshifts and scale factors if snapinput is set.
    
    Parameters
    ----------
    snapinput : int or array, optional
        Subset of snapnums for which redshifts and scale factors are desired
    
    Returns
    -------
    (snaps,) redshifts, times : ndarrays
        Arrays of redshifts and scale factors corresponding to snapnums
    """

    snapfile = pkg_resources.resource_filename('indratools', 'data/FFT_snapinfo.txt')
    snaps,redshifts,times = np.loadtxt(snapfile,unpack=True)

    if snapinput == None:
        return snaps, redshifts, times
    else:
        return redshifts[snapinput],times[snapinput]


def get_pklin():
    """Reads and returns the k and P(k) values of the CAMB linear power spectrum, 
    normalized to z=0, for the Indra cosmology.
    
    Returns
    -------
    k, pk : ndarrays
        Arrays of wave-vector and power spectrum values
    """

    pkfile = pkg_resources.resource_filename('indratools','data/pk_indra7313_CAMB.txt')
    k,pk = np.loadtxt(pkfile,unpack=True)
    
    return k, pk

