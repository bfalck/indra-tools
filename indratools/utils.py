"""
Utility functions for the Indra suite of simulations hosted on the SciServer.

Written by Bridget Falck, 2018-2020



Inputs: 
- ``snapinput`` is an optional subset of snapshots for which redshifts and scale factors
    are desired for particle or FFT snapshots.
- ``runid`` specifies the Indra run, and is ignored if a ``datadir`` is specified; 
    it is either an integer from 0 to 511 or a tuple containing (``X``,``Y``,``Z``)



Classes
-------

Run(runid)
    Define the runnum (if runid is a tuple of (X,Y,Z)) and X, Y, Z values (if runid is an int)
    that specify this Indra run.

Methods
-------

get_run_num(x,y,z)
    Helper function to return the raveled run number (0 to 511) from the unraveled
    x, y, z identifiers (each 0 to 7).
get_xyz(run_num)
    Helper function to return the x, y, z identifiers from the run number (0 to 511).
fdb_snaps(runnum=None,Container=True,Dask=False,DaskLocal=False)
    Returns snapshots on FileDB: full sets (default), or those in directory of runnum.
    Optional boolean inputs specify the environment to determine filesystem paths.
fdb_runs(loc=None,getstring=True,Container=True,Dask=False,DaskLocal=False)
    Returns Indra runs (as string 'X_Y_Z' (default) or as a runnum) on FileDB: full
    runs (default) or runs located in FileDB node specified by loc = node_vol, node_num.
    For example, filedb/data05_01 specified by loc = (5,1).
    Optional boolean inputs specify the environment to determine filesystem paths.
fdb_loc(runid)
    Returns file location of given Indra run on the FileDB data volumes, as 
    mounted on SciServer Compute containers.
get_loc(runid,snapnum=None)
    Determines directory from which to read simulation snapshot file.
    Defaults to FileDB location for all halo and FFT data, and for all snapshots
    stored on FileDB. Otherwise it returns the datascope location based on ``ds_basedir``.
part_snapinfo(snapinput=None)
    Reads and returns all snapnums, redshifts, and scale factors for all 64 particle 
    snapshots, or a subset of redshifts and scale factors if snapinput is set.
FFT_snapinfo(snapinput=None)
    Reads and returns all snapnums, redshifts, and scale factors for all 505 FFT 
    snapshots, or a subset of redshifts and scale factors if snapinput is set.
get_pklin()
    Reads and returns the k and P(k) values of the CAMB linear power spectrum, 
    normalized to z=0, for the Indra cosmology.
cic_pk(snapnum,ngrid=512,nruns=320)
    Reads and returns the power spectra measured from CIC interpolation of
    the particle positions of multiple runs and one snapshot.

"""

import numpy as np
import pkg_resources
import glob


ds_basedir = '/home/idies/workspace/indra_dss/'


def get_run_num(x,y,z):
    '''Helper function to figure out raveled index from unraveled index.

    Parameters
    ----------
    x, y, z : int
        The integers specifying run x_y_z, each from 0 to 7

    Returns
    -------
    int
        The ID of the run as an integer, from 0 to 511
    '''
#    return x*64+y*8+z
    return np.ravel_multi_index((x,y,z),(8,8,8))


def get_xyz(run_num):
    '''Helper function to figure out unraveled index from raveled index.

    Parameters
    ----------
    run_num : int
        The ID of the run as an integer, from 0 to 511
    
    Returns
    -------
    x, y, z : ints
        Integers from 0 to 7 identifying run x_y_z
    '''
#    return run_num//64, run_num//8 % 8, run_num % 8
    return np.unravel_index(run_num,(8,8,8))


class Run:
    """
    Specifies the current Indra simulation as both a number (from 0 to 511) 
    and 3 integers, X_Y_Z (each go from 0 to 7), corresponding to the 
    raveled and unraveled indices of an 8x8x8 cube. Instantiated with either
    the number or the 3 integers as a tuple.
    
    Attributes
    ----------
    num : int
        The ID of the run as an integer, from 0 to 511
    X : int
        The first integer of run X_Y_Z, from 0 to 7
    Y : int
        The second integer of run X_Y_Z, from 0 to 7
    Z : int
        The third integer of run X_Y_Z, from 0 to 7

    """

    def __init__(self, runid):
        """
        Parameters
        ----------
        runid : int or tuple
            Specifies the Indra run either as an integer from 0 to 511
            or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
            where X, Y, and Z each go from 0 to 7.
        """
        if isinstance(runid, int) or isinstance(runid, np.integer):
            self.num = runid
            self.X, self.Y, self.Z = get_xyz(runid)
        elif isinstance(runid, tuple):
            self.X, self.Y, self.Z = runid
            self.num = get_run_num(self.X,self.Y,self.Z)

    
def fdb_snaps(runid=None,Container=True,Dask=False,DaskLocal=False):
    """Return list of snapshots on FileDB by running glob on the filesystem, either
    for a specific run given by runnum, or return list of the full sets of snapshots
    that are on FileDB (for every run).
    
    Parameters
    ----------
    runid : int or tuple (default None)
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    Container : boolean (default True)
        Check filesystem as mounted on a Compute container.
    Dask : boolean (default False)
        Check filesystem as mounted on the Dask cluster.
    DaskLocal : boolean (default False)
        Check filesystem as seen by Dask worker on local FileDB node.

    Returns
    -------
    ndarray
        Sorted array of snapshots as integers.
    """

    if (Dask or DaskLocal):
        Container = False

    if Container:
        basedir='/home/idies/workspace/indra_filedb/'
    
        
    if runid is None:
        # make sure to pick a non-full run: choosing 2_0_1
        if Container:
            datadir = basedir+'data03_01/2_0_1/'
        elif Dask:
            datadir = '/sciserver/filedb03-01/cosmo/indra/2_0_1/'
        elif DaskLocal:
            datadir = '/srv/data01/cosmo/indra/2_0_1/'

        globstr = datadir+'snapdir_*'
    else:
        # first get X_Y_Z from runnum, then use it to list snaps; filedb node can be wildcard
        run = Run(runid)
        x, y, z = (run.X,run.Y,run.Z)
        runstr = f'{x}_{y}_{z}'
        if Container:
            globstr = basedir+'*/'+runstr+'/snapdir_*'
        elif Dask:
            globstr = '/sciserver/*/cosmo/indra/'+runstr+'/snapdir_*'
        elif DaskLocal:
            globstr = '/srv/*/cosmo/indra/'+runstr+'/snapdir_*'

    snaps = sorted([np.int(line[-2:]) for line in glob.iglob(globstr)])
    return np.array(snaps)
            

def fdb_runs(loc=None,getstring=True,Container=True,Dask=False,DaskLocal=False):
    """Returns list of full runs (containing all snapshots) on FileDB.
    If loc is set, return list of runs at that FileDB location (one of 36 volumes).
    
    Parameters
    ----------
    loc : (node_volume, node_number) (optional)
        Tuple that identifies one of the 36 FileDB locations, where 
        node_volume is an int from 1 to 12, node_number an int from 1 to 3.
        If set, return list of runs at that location. If not, return list
        of full runs (containing every snapshot) on FileDB.
    getstring : boolean (default True)
        If True, return runs as string identifiers ('X_Y_Z').
        If False, return as list of integers.
    Container : boolean (default True)
        Check filesystem as mounted on a Compute container.
    Dask : boolean (default False)
        Check filesystem as mounted on the Dask cluster.
    DaskLocal : boolean (default False)
        Check filesystem as seen by Dask worker on local FileDB node.

    
    Returns
    -------
    list
        Sorted list of string run identifiers ('X_Y_Z') or integer run numbers.
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
    

def fdb_loc(runid):    
    """Helper function to find file location of given Indra run on the FileDB
    data volumes, as mounted on SciServer Compute containers.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    
    Returns
    -------
    filename : string
        Full path of the directory containing the run specified by runid.
    """

    run = Run(runid)

    # Get list of filedb locations: start with 08-01
    fd = []
    for f in range(8,13):
        for d in range(1,4):
            fd.append('/home/idies/workspace/indra_filedb/data{:02d}_{:02d}/'.format(f,d))
    for f in range(1,8):
        for d in range(1,4):
            fd.append('/home/idies/workspace/indra_filedb/data{:02d}_{:02d}/'.format(f,d))

    return fd[run.num % 36]+'{}_{}_{}/'.format(run.X,run.Y,run.Z)


def get_loc(runid,snapnum=None):
    """Determines directory from which to read simulation snapshot file as mounted on
    SciServer compute containers.
    
    If snapnum is not provided, returns FileDB location of input run. Otherwise, checks
    whether provided runid and snapnum are on FileDB. If not, returns the DataScope 
    location based on ``ds_basedir``, a global parameter here and in read_indra.py.
    Thus it returns the FileDB location for all halo catalogs and FFT data and for 
    all snapshots that are stored on FileDB.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int (optional)
        Only required to determine location of snapshot files.

    Returns
    -------
    filename : string
        Full path of the directory containing the run specified by runid.
    """

    if snapnum is None:
        return fdb_loc(runid)
    else:
        run = Run(runid)
        if run.num in fdb_runs(getstring=False) or snapnum in fdb_snaps(runid):
            return fdb_loc(runid)
        else:
            return f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'


def part_snapinfo(snapinput=None):
    """Reads and returns all snapnums, redshifts, and scale factors for all 64 particle
    snapshots, or a subset of redshifts and scale factors if snapinput is set.
    
    Parameters
    ----------
    snapinput : int or array_like, optional
        Subset of snapnums for which redshifts and scale factors are desired
    
    Returns
    -------
    (snaps,) redshifts, times : ndarrays
        Arrays of redshifts and scale factors corresponding to snapnums
    """
    
    snapfile = pkg_resources.resource_filename('indratools', 'data/snapinfo.txt')
    snaps,redshifts,times = np.loadtxt(snapfile,unpack=True)
    
    if snapinput is None:
        return snaps, redshifts, times
    else:
        return redshifts[snapinput],times[snapinput]
    
    
def FFT_snapinfo(snapinput=None):
    """Reads and returns all snapnums, redshifts, and scale factors for all 505 FFT
    snapshots, or a subset of redshifts and scale factors if snapinput is set.
    
    Parameters
    ----------
    snapinput : int or array_like, optional
        Subset of snapnums for which redshifts and scale factors are desired
    
    Returns
    -------
    (snaps,) redshifts, times : ndarrays
        Arrays of redshifts and scale factors corresponding to snapnums
    """

    snapfile = pkg_resources.resource_filename('indratools', 'data/FFT_snapinfo.txt')
    snaps,redshifts,times = np.loadtxt(snapfile,unpack=True)

    if snapinput is None:
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


def cic_pk(snapnum,ngrid=512,nruns=320):
    """Reads and returns the power spectra measured from CIC interpolation of
    the particle positions of multiple runs and one snapshot.
    
    Parameters
    ----------
    snapnum : int
        snapshot number
    ngrid : int (default 512)
        number of grid cells used for CIC interpolation (determines Nyquist frequency)
    nruns : int (default 320)
        number of Indra runs included (e.g. 320 is for runs 128 to 447)
    
    Returns
    -------
    dict of {'k': ndarray, 'ps' : ndarray, 'poisson' : ndarray}
        'k' is a 1-d array of bin values, 'ps' is a [len(k),nruns] array of P(k) values,
        and 'poisson' is a 1-d array of 1/sqrt(N) for N counts in each bin.
    """
    
    basedir = '/home/idies/workspace/indra/pk/' 
    k = np.load(f'{basedir}/psbins_ng{ngrid}.npy')
    pk = np.load(f'{basedir}/psvals_ng{ngrid}_s{snapnum}_{nruns}runs.npy')
    err = np.load(f'{basedir}/pserr_ng{ngrid}.npy')
    
    return {'k':k, 'ps': pk, 'poisson': err}
    
    