"""
Utility functions for the Indra suite of simulations hosted on the SciServer.

Written by Bridget Falck, 2018-2020


TO DO: documentation for fdb functions


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
fdb_snaps(runnum=None,...)
    Returns snapshots on FileDB: full sets (default), or those in directory of runnum
fdb_runs(loc=None,getstring=True,...)
    Returns Indra runs (as string 'X_Y_Z' (default) or runnum) on FileDB: full runs
    (default) or runs located in FileDB node specified by loc = node_vol, node_num.
    E.g., filedb/data05_01 specified by loc = (5,1).

"""

import numpy as np
import pkg_resources
import glob


def get_run_num(x,y,z):
    '''Helper function to figure out raveled index from unraveled index'''
#    return x*64+y*8+z
    return np.ravel_multi_index((x,y,z),(8,8,8))

def get_xyz(run_num):
    '''Helper function to figure out unraveled index from raveled index'''
#    return run_num//64, run_num//8 % 8, run_num % 8
    return np.unravel_index(run_num,(8,8,8))


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


    
def fdb_snaps(runnum=None,Container=True,Dask=False,DaskLocal=False):
    """Return list of snapshots on FileDB.
    If no runnum given, pick a non-full run to choose snaps that exist for every run
    """
    if (Dask or DaskLocal):
        Container = False

    if Container:
        basedir='/home/idies/workspace/indra_filedb/'

    if runnum is None:
        # make sure to pick a non-full run: choosing 2_0_1
        if Container:
            datadir = basedir+'data03_01/2_0_1/'
        elif Dask:
            datadir = '/sciserver/filedb03-01/cosmo/indra/2_0_1/'
        elif DaskLocal:
            datadir = '/srv/data01/cosmo/indra/2_0_1/'

        globstr = datadir+'snapdir_*'
    else:
        # first get runid from runnum, then use it to list snaps; filedb node can be wildcard
        x, y, z = get_xyz(runnum)
        runid = f'{x}_{y}_{z}'
        if Container:
            globstr = basedir+'*/'+runid+'/snapdir_*'
        elif Dask:
            globstr = '/sciserver/*/cosmo/indra/'+runid+'/snapdir_*'
        elif DaskLocal:
            globstr = '/srv/*/cosmo/indra/'+runid+'/snapdir_*'

    snaps = sorted([np.int(line[-2:]) for line in glob.iglob(globstr)])
    return np.array(snaps)
            

def fdb_runs(loc=None,getstring=True,Container=True,Dask=False,DaskLocal=False):
    """Returns list of full runs (all snaps) on FileDB.
    If some input set, return list of runs at that loc
    INPUT: loc is a tuple of (node volume, node num)
    """
    if (Dask or DaskLocal):
        Container = False

    if Container:
        basedir='/home/idies/workspace/indra_filedb/'

    if loc is None:
        if Container:
            globstr = basedir+'*/*/snapdir_002'
        elif Dask:
            globstr = '/sciserver/*/cosmo/indra/*/snapdir_002'
        elif DaskLocal:
            globstr = '/srv/*/cosmo/indra/*/snapdir_002'

        runs = sorted([line[-17:-12] for line in glob.iglob(globstr)])

    else:
        node_vol, node_num = loc
        if Container:
            datadir = f"{basedir}data{node_vol:02}_{node_num:02}/"
        elif Dask:
            datadir = f"/sciserver/filedb{node_vol:02}-{node_num:02}/cosmo/indra/"
        elif DaskLocal:
            datadir = f"/srv/data{node_num:02}/cosmo/indra/"

        runs = sorted([line[-5:] for line in glob.iglob(datadir+'*')])

    run_nums = [get_run_num(np.int(runstr[0]),np.int(runstr[2]),np.int(runstr[4])) for runstr in runs]
    if getstring:
        return runs
    else:
        return run_nums
    
    
    
    