"""
Reading functions for the Indra suite of simulations hosted on SciServer.

**Snapshot functions:**

- getheader: Returns a dictionary of header parameters.
- getpos: Returns an array of particle positions of shape ``[1024**3,3]``
    from one full snapshot.
- getparticles: Returns particle positions, velocities, and IDs: ``pos, vel, ids = getparticles(...)``.
    `pos` and `vel` are arrays of shape ``[1024**3,3]``, and `ids` have shape ``[1024**3]``.

**Halo and subhalo functions:**

- getfofheader: Reads the header of the FOF data. Returns the total number of FOF 
    groups contained in all files of given snapshot.
- getsubheader: Reads all headers of the SUBFIND files and returns the total number
    of subhalos (default) or FOF groups (if `getfof`=True) for the given snapshot.
- getfof: Reads the number of particles in each FOF group. If `getOffset`=True, also returns
    the Offset array needed to index the IDs of the member particles: 
    ``groupLen, groupOffset = getfof(...)``
- getfofids: Reads the `groupLen`, `groupOffset`, and particle ID arrays for the FOF groups.
- getsubcat: Reads the SUBFIND subhalo catalogs and returns a dictionary of values. Some halo 
    properties are defined for each FOF group, and some for each subhalo. The
    dictionary contains the `subLen` and `subOffset` arrays needed to index the IDs
    of the subhalo member particles (as for the FOF groups).
- getsubids: Reads and returns the IDs of the particles in each subhalo.

**FFT functions:**

- getfft: Reads and returns the real and imaginary components of the Fourier-space density
    grid at output `tnum` and returns the scalefactor of this `tnum`.
- getkvals: Builds and returns the *x*, *y*, and *z* components of the *k*-vectors associated
    with the FFT data. Each are arrays with the same shapes as `fft_re` and `fft_im`.
"""

# package: from .utils import *
# dev: from utils import *
from .utils import *
import numpy as np

ds_basedir = '/home/idies/workspace/indra_dss/'

# Run/snap combinations with missing subfiles where there are no subhalos anyway:
MissingFiles_OK = [(3,0,1,1)]
# Run/snap combinations with missing/empty subfiles where there should be some:
MissingFiles_Problem = []




def _readheader(f):
    """Utility function for reading Indra snapshot files"""

    # first sanity check: header_size
    if np.fromfile(f,np.int32,1) != 256:
        return None
    header={}
        
    npall = np.fromfile(f,np.int32,6) # DM only simulations - other np's and masses are 0
    header['np_file'] = npall[1] # number of particles in this file (needed to read positions)
    massall = np.fromfile(f,np.float64,6)
    header['mass'] = massall[1]*1.e10 # now in units of Msun/h
    header['time'],header['redshift'] = np.fromfile(f,np.float64,2)
    header['flag_sfr'],header['flag_feedback'] = np.fromfile(f,np.int32,2)
    nptotal = np.fromfile(f,np.int32,6)
    header['npart'] = nptotal[1]
    header['flag_cooling'],header['num_files'] = np.fromfile(f,np.int32,2)
    header['BoxSize'],header['omega0'],header['omegaL'],header['hubble'] = np.fromfile(f,np.float64,4)
    header['flag_stellarage'],header['flag_metals'],header['hashtabsize'] = np.fromfile(f,np.int32,3)
    # read rest of header_size + 2 filler integers (so next read can be particle positions):
    empty = np.fromfile(f,np.int32,23)
    return header


def getheader(runid=128,snapnum=0,filenum=0,datadir=None,datascope=False,verbose=False):
    """Reads the Gadget header of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 128 to 511
        or as a length 3 tuple giving the 3-digit ID as (`X`,`Y`,`Z`)
        where `X` goes from 2 to 7 and `Y` and `Z` each go from 0 to 7.
        Ignored if `datadir` is set.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation `X_Y_Z`.
        Default is to read from the output of ``get_loc(runid,snapnum)``.
    datascope : boolean, optional
        If True, read from ``/datascope_path/indraX/X_Y_Z/`` (default False).
        Ignored if `datadir` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    header : dict
        Dictionary object containing header information.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid,snapnum)

    if (verbose == True): print('Reading from {}'.format(datadir))

    filename = '{0}/snapdir_{1:03d}/snapshot_{1:03d}.'.format(datadir,snapnum)


    try:
        with open(filename+str(filenum),'rb') as f:
            header = _readheader(f)
    except FileNotFoundError:
        # doesn't exist
        print('No file {}: returning None'.format(filename+str(0)))
        return None

    return header

def _readpos(f,npfile):
    """Utility function for `getpos()`"""
    
    thispos = np.fromfile(f,np.float32,3*npfile)
    thispos = np.reshape(thispos, [npfile, 3])
    
    return thispos

def _readsnap(f,npfile):
    """Utility function for `getparticles()`"""
    
    thispos = _readpos(f,npfile)
    empty = np.fromfile(f,np.int32,2)
    thisvel = np.fromfile(f,np.float32,3*npfile)
    thisvel = np.reshape(thisvel, [npfile, 3])
    empty = np.fromfile(f,np.int32,2)
    thisID = np.fromfile(f,np.int64,npfile)
    
    return thispos,thisvel,thisID

def getpos(runid,snapnum,datadir=None,datascope=False,verbose=False):
    """Reads the particle positions (only) of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 128 to 511
        or as a length 3 tuple giving the 3-digit ID as (`X`,`Y`,`Z`)
        where `X` goes from 2 to 7 and `Y` and `Z` each go from 0 to 7.
        Ignored if `datadir` is set.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation `X_Y_Z`.
        Default is to read from the output of ``get_loc(runid,snapnum)``.
    datascope : boolean, optional
        If True, read from ``/datascope_path/indraX/X_Y_Z/`` (default False).
        Ignored if `datadir` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    pos : ndarray
        Particle positions (in Mpc/*h*) in array of shape ``[1024**3,3]``
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid,snapnum)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    filename = '{0}/snapdir_{1:03d}/snapshot_{1:03d}.'.format(datadir,snapnum)

    # Check that files exist:
    try:
        with open(filename+str(0),'rb') as f:
            pass
    except FileNotFoundError:
        # doesn't exist
        print('No file {}: returning None'.format(filename+str(0)))
        return None
    
    # indra7: for some runs and snapshots, NTask is different
    # Need to read from header['num_files']
    header = getheader(runid,snapnum,datadir=datadir,datascope=datascope)
    NTask = header['num_files']
    nparticles = 1024
    
    # loop through files
    pos = np.empty((nparticles**3,3),np.float32)
    istart = 0
    for i in np.arange(0,NTask):
        f = open(filename+str(i),'rb')
        
        header = _readheader(f)
        thispos = _readpos(f,header['np_file'])
        
        f.close()
    
        pos[istart:(istart+header['np_file']),:] = thispos
        istart = istart + header['np_file']

    return pos

def getparticles(runid,snapnum,datadir=None,datascope=False,sort=False,verbose=False):
    """Reads the particle positions, velocities, and IDs of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 128 to 511
        or as a length 3 tuple giving the 3-digit ID as (`X`,`Y`,`Z`)
        where `X` goes from 2 to 7 and `Y` and `Z` each go from 0 to 7.
        Ignored if `datadir` is set.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation `X_Y_Z`.
        Default is to read from the output of ``get_loc(runid,snapnum)``.
    datascope : boolean, optional
        If True, read from ``/datascope_path/indraX/X_Y_Z/`` (default False).
        Ignored if `datadir` is set.
    sort : boolean, optional
        If True, sort the positions and velocities by the particle IDs.
        The sorting is done on-the-fly to minimize reading time.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    pos, vel, ids : tuple of ndarrays
        Particle positions (in Mpc/*h*) and velocities (in km/s) in arrays of 
        shape ``[1024**3,3]`` and IDs in array of shape ``[1024**3]``. If `sort`=True,
        ID array is equivalent to ``numpy.arange(1024**3)``.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid,snapnum)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    filename = '{0}/snapdir_{1:03d}/snapshot_{1:03d}.'.format(datadir,snapnum)

    # Check that files exist:
    try:
        with open(filename+str(0),'rb') as f:
            pass
    except FileNotFoundError:
        # doesn't exist
        print('No file {}: returning None'.format(filename+str(0)))
        return None,None,None

    #NTask = 256
    # indra7: for some runs and snapshots, NTask is different
    header = getheader(runid,snapnum,datadir=datadir,datascope=datascope)
    NTask = header['num_files']
    nparticles = 1024

    # loop through files
    pos = np.empty((nparticles**3,3),np.float32)
    vel = np.empty((nparticles**3,3),np.float32)
    ids = np.empty((nparticles**3),np.int64)
    istart = 0
    for i in np.arange(0,NTask):
        f = open(filename+str(i),'rb')
        header = _readheader(f)
        thispos,thisvel,thisID = _readsnap(f,header['np_file'])
        f.close()

        thisID -= 1 # make IDs start at 0 to use to index pos and vel
        
        if sort:
            pos[thisID] = thispos
            vel[thisID] = thisvel*np.sqrt(header['time'])
        else:
            pos[istart:(istart+header['np_file']),:] = thispos
            vel[istart:(istart+header['np_file']),:] = thisvel*np.sqrt(header['time'])
            ids[istart:(istart+header['np_file'])] = thisID
    
        istart = istart + header['np_file']
        
    if sort: ids = np.arange(header['npart'])

    return pos,vel,ids


"""
Halo and subhalo functions
"""

def getfofheader(runid,snapnum,datadir=None,datascope=False,verbose=False):
    """Reads the FOF header (total number of groups) of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 128 to 511
        or as a length 3 tuple giving the 3-digit ID as (`X`,`Y`,`Z`)
        where `X` goes from 2 to 7 and `Y` and `Z` each go from 0 to 7.
        Ignored if `datadir` is set.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation `X_Y_Z`.
        Default is to read from the output of ``get_loc(runid,snapnum)``.
    datascope : boolean, optional
        If True, read from ``/datascope_path/indraX/X_Y_Z/`` (default False).
        Ignored if `datadir` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    TotNgroups : int
        Total number of FOF groups in this snapshot.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid,snapnum)

    if (verbose == True): print('Reading from {}'.format(datadir))
        
    tabfile = '{0}/snapdir_{1:03d}/group_tab_{1:03d}.'.format(datadir,snapnum)

    try:
        with open(tabfile+str(0),'rb') as f:
            Ngroups, Nids, TotNgroups, NTask = np.fromfile(f, np.int32, 4)
    except FileNotFoundError:
        # doesn't exist
        print('No file {}'.format(tabfile+str(0)))
        TotNgroups = 0

    return TotNgroups

def getfof(runid,snapnum,datadir=None,datascope=False,getOffset=False,verbose=False):
    """Reads the FOF group info of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 128 to 511
        or as a length 3 tuple giving the 3-digit ID as (`X`,`Y`,`Z`)
        where `X` goes from 2 to 7 and `Y` and `Z` each go from 0 to 7.
        Ignored if `datadir` is set.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation `X_Y_Z`.
        Default is to read from the output of ``get_loc(runid,snapnum)``.
    datascope : boolean, optional
        If True, read from ``/datascope_path/indraX/X_Y_Z/`` (default False).
        Ignored if `datadir` is set.
    getOffset : boolean, optional
        If True, returns `groupLen`, `groupOffset`, otherwise just `groupLen` (default False)
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    groupLen (, groupOffset) : ndarray or tuple of ndarrays
        Number of particles in each FOF group, and if `getOffset` = True,
        Offset array needed in order to index the IDs of the group's 
        constituent particles.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid,snapnum)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    tabfile = '{0}/snapdir_{1:03d}/group_tab_{1:03d}.'.format(datadir,snapnum)

    # Need NTask (number of files) in case it's not 256, and TotNgroups
    try:
        with open(tabfile+str(0),'rb') as f:
            Ngroups, Nids, TotNgroups, NTask = np.fromfile(f, np.int32, 4)
    except FileNotFoundError:
        # doesn't exist
        print('No file {}'.format(tabfile+str(0)))
        NTask = 0
        TotNgroups = 0

    # Don't read if no groups...
    if TotNgroups == 0: 
        print("No FOF groups in {}_{}_{} snapshot {}: returning None".format(run.X,run.Y,run.Z,snapnum))
        if getOffset:
            return None,None
        else:
            return None
    elif getOffset:
        groupLen = np.zeros(TotNgroups,dtype=np.int32)
        groupOffset = np.zeros(TotNgroups,dtype=np.int32)

        istartGroup = 0
        istartIDs = 0
        for i in np.arange(0,NTask):
            f = open(tabfile+str(i), 'rb')
    
            Ngroups, Nids, TotNgroups, NTask = np.fromfile(f, np.int32, 4)
            if Ngroups > 0:
                locLen = np.fromfile(f,np.int32,Ngroups)
                locOffset = np.fromfile(f,np.int32,Ngroups)
                groupLen[istartGroup:(istartGroup+Ngroups)] = locLen
                groupOffset[istartGroup:(istartGroup+Ngroups)] = locOffset+istartIDs
                istartGroup += Ngroups
                istartIDs += Nids
            f.close()
    else:
        groupLen = np.zeros(TotNgroups,dtype=np.int32)

        istartGroup = 0
        for i in np.arange(0,NTask):
            f = open(tabfile+str(i), 'rb')
    
            Ngroups, Nids, TotNgroups, NTask = np.fromfile(f, np.int32, 4)
            if Ngroups > 0:
                locLen = np.fromfile(f,np.int32,Ngroups)
                groupLen[istartGroup:(istartGroup+Ngroups)] = locLen
                istartGroup += Ngroups
            f.close()

    if getOffset:
        return groupLen,groupOffset
    else:
        return groupLen

def getfofids(runid,snapnum,datadir=None,datascope=False,verbose=False):
    """Reads the FOF group info (including particle IDs) of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 128 to 511
        or as a length 3 tuple giving the 3-digit ID as (`X`,`Y`,`Z`)
        where `X` goes from 2 to 7 and `Y` and `Z` each go from 0 to 7.
        Ignored if `datadir` is set.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation `X_Y_Z`.
        Default is to read from the output of ``get_loc(runid,snapnum)``.
    datascope : boolean, optional
        If True, read from ``/datascope_path/indraX/X_Y_Z/`` (default False).
        Ignored if `datadir` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    groupLen, groupOffset, groupIDs : tuple of ndarrays
        Number of particles in each FOF group, Offset array needed
        in order to index the ID array, and the IDs of each group's
        constituent particles.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid,snapnum)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    idsfile = '{0}/snapdir_{1:03d}/group_ids_{1:03d}.'.format(datadir,snapnum)
 
    # loop through NTask files
    try:
        with open(idsfile+str(0),'rb') as f:
            Ngroups, Nids, TotNgroups, NTask = np.fromfile(f, np.int32, 4)
    except FileNotFoundError:
        # doesn't exist
        print('No file {}'.format(tabfile+str(0)))
        NTask = 0
        TotNgroups = 0

    # Don't read if no groups...
    if TotNgroups == 0: 
        print("No FOF groups in {}_{}_{} snapshot {}: returning None".format(run.X,run.Y,run.Z,snapnum))
        return None,None,None
    else:
        groupLen,groupOffset = getfof(runid,snapnum,datadir,getOffset=True)
        TotNids = np.sum(groupLen,dtype=np.int64)
        groupids = np.zeros(TotNids,dtype=np.int64)

        istart = 0
        for i in np.arange(0,NTask):
            f = open(idsfile+str(i),'rb')

            Ngroups, Nids, TotNgroups, NTask = np.fromfile(f, np.int32, 4)
            if Nids > 0:
                locIDs = np.fromfile(f,np.int64,Nids)
                # bitshift to remove the hash table info from the IDs
                groupids[istart:(istart+Nids)] = np.bitwise_and(locIDs[:], (np.int64(1)<<34) - 1)
                istart += Nids
#            else: print('No groups in file %d' % i)
            f.close()
    
    return groupLen,groupOffset,groupids-1

def getsubheader(runid,snapnum,datadir=None,datascope=False,getfof=False,verbose=False):
    """Reads the SUBFIND header (total number of subhalos) of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 128 to 511
        or as a length 3 tuple giving the 3-digit ID as (`X`,`Y`,`Z`)
        where `X` goes from 2 to 7 and `Y` and `Z` each go from 0 to 7.
        Ignored if `datadir` is set.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation `X_Y_Z`.
        Default is to read from the output of ``get_loc(runid,snapnum)``.
    datascope : boolean, optional
        If True, read from ``/datascope_path/indraX/X_Y_Z/`` (default False).
        Ignored if `datadir` is set.
    getfof : boolean, optional 
        If true, return `TotNgroups` instead of `TotNsubs` (default False).
        (For example, when group_tab files are not available).
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    TotNsubs (or TotNgroups if `getfof` is True) : int
        Total number of subhalos or FOF groups in this snapshot.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid) #snapnum not specified, all files are on FileDB

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    tabfile = '{0}/postproc_{1:03d}/sub_tab_{1:03d}.'.format(datadir,snapnum)

    # don't read if there are missing files for any reason:
    # (because MissingFiles_OK is at early snapshot where there are no halos anyway)
    if (run.X,run.Y,run.Z,snapnum) in MissingFiles_OK+MissingFiles_Problem:
        return 0

    # open first file to get NTask just in case (but it should be 256) and TotNgroups
    f = open(tabfile+str(0),'rb')
    Ngroups,Nids,TotNgroups,NTask = np.fromfile(f,np.int32,4)
    f.close()

    if getfof: return TotNgroups
    else:
        TotNsubs = 0
        for i in np.arange(0,NTask):
            f = open(tabfile+str(i),'rb')
    #        print('opening file '+tabfile+str(i))
            Ngroups, Nids, TotNgroups, NTask, Nsubs = np.fromfile(f, np.int32, 5)
            TotNsubs += Nsubs
            f.close()

        return TotNsubs

def getsubcat(runid,snapnum,datadir=None,datascope=False,verbose=False):
    """Reads the subhalo catalog of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 128 to 511
        or as a length 3 tuple giving the 3-digit ID as (`X`,`Y`,`Z`)
        where `X` goes from 2 to 7 and `Y` and `Z` each go from 0 to 7.
        Ignored if `datadir` is set.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation `X_Y_Z`.
        Default is to read from the output of ``get_loc(runid,snapnum)``.
    datascope : boolean, optional
        If True, read from ``/datascope_path/indraX/X_Y_Z/`` (default False).
        Ignored if `datadir` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    catalog : dictionary
        Dictionary of halo and subhalo properties and indexing information.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid) #snapnum not specified, all files are on FileDB

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    tabfile = '{0}/postproc_{1:03d}/sub_tab_{1:03d}.'.format(datadir,snapnum)
    
    # open first file to get NTask just in case (but it should be 256) and TotNgroups
    f = open(tabfile+str(0),'rb')
    Ngroups,Nids,TotNgroups,NTask = np.fromfile(f,np.int32,4)
    f.close()
    TotNsubs = getsubheader(runid,snapnum,datadir)

#    if TotNgroups == 0: return None
    if TotNsubs == 0: 
        if (run.X,run.Y,run.Z,snapnum) in MissingFiles_Problem:
            print("Missing data in {}_{}_{} snapshot {}: returning None".format(run.X,run.Y,run.Z,snapnum))
        else:
            print("No subhalos in {}_{}_{} snapshot {}: returning None".format(run.X,run.Y,run.Z,snapnum))

        return None

    else:
        # Initialize data containers
        catalog = {}
        catalog['NsubPerHalo'] = np.zeros(TotNgroups,dtype=np.int32)
        catalog['FirstSubOfHalo'] = np.zeros(TotNgroups,dtype=np.int32) # file specific!
    
        catalog['subLen'] = np.zeros(TotNsubs,dtype=np.int32)
        catalog['subOffset'] = np.zeros(TotNsubs,dtype=np.int32) # file specific!
        catalog['subParentHalo'] = np.zeros(TotNsubs,dtype=np.int32) # file specific!
    
        catalog['M200mean'] = np.zeros(TotNgroups,dtype=np.float32)
        catalog['R200mean'] = np.zeros(TotNgroups,dtype=np.float32)
        catalog['M200crit'] = np.zeros(TotNgroups,dtype=np.float32)
        catalog['R200crit'] = np.zeros(TotNgroups,dtype=np.float32)
        catalog['M200tophat'] = np.zeros(TotNgroups,dtype=np.float32)
        catalog['R200tophat'] = np.zeros(TotNgroups,dtype=np.float32)
    
        catalog['SubPos'] = np.zeros((TotNsubs,3),dtype=np.float32)
        catalog['SubVel'] = np.zeros((TotNsubs,3),dtype=np.float32)
        catalog['SubVelDisp'] = np.zeros(TotNsubs,dtype=np.float32)
        catalog['SubVmax'] = np.zeros(TotNsubs,dtype=np.float32)
        catalog['SubSpin'] = np.zeros((TotNsubs,3),dtype=np.float32)
        catalog['SubMostBoundID'] = np.zeros(TotNsubs,dtype=np.int64)
        catalog['SubHalfMass'] = np.zeros(TotNsubs,dtype=np.float32)
        
        # loop over numfiles (NTask)
        istartSub = 0
        istartGroup = 0
        istartIDs = 0
        for i in np.arange(0,NTask):
            f = open(tabfile+str(i),'rb')
            Ngroups, Nids, TotNgroups, NTask, Nsubs = np.fromfile(f, np.int32, 5)
            
            if Nsubs > 0:
            # Read catalog: indexes need to include offset from previous files (istarts)
                catalog['NsubPerHalo'][istartGroup:(istartGroup+Ngroups)] = np.fromfile(f,np.int32,Ngroups)
                catalog['FirstSubOfHalo'][istartGroup:(istartGroup+Ngroups)] = np.fromfile(f,np.int32,Ngroups)+istartSub
       
                catalog['subLen'][istartSub:(istartSub+Nsubs)] = np.fromfile(f,np.int32,Nsubs)
                catalog['subOffset'][istartSub:(istartSub+Nsubs)] = np.fromfile(f,np.int32,Nsubs)+istartIDs
                catalog['subParentHalo'][istartSub:(istartSub+Nsubs)] = np.fromfile(f,np.int32,Nsubs)+istartGroup
                
                catalog['M200mean'][istartGroup:(istartGroup+Ngroups)] = np.fromfile(f,np.float32,Ngroups)*1.e10
                catalog['R200mean'][istartGroup:(istartGroup+Ngroups)] = np.fromfile(f,np.float32,Ngroups)
                catalog['M200crit'][istartGroup:(istartGroup+Ngroups)] = np.fromfile(f,np.float32,Ngroups)*1.e10
                catalog['R200crit'][istartGroup:(istartGroup+Ngroups)] = np.fromfile(f,np.float32,Ngroups)
                catalog['M200tophat'][istartGroup:(istartGroup+Ngroups)] = np.fromfile(f,np.float32,Ngroups)*1.e10
                catalog['R200tophat'][istartGroup:(istartGroup+Ngroups)] = np.fromfile(f,np.float32,Ngroups)

                thisxyz = np.fromfile(f,np.float32,3*Nsubs)
                catalog['SubPos'][istartSub:(istartSub+Nsubs),:] = np.reshape(thisxyz,[Nsubs,3])
                thisxyz = np.fromfile(f,np.float32,3*Nsubs)
                catalog['SubVel'][istartSub:(istartSub+Nsubs),:] = np.reshape(thisxyz,[Nsubs,3])
                catalog['SubVelDisp'][istartSub:(istartSub+Nsubs)] = np.fromfile(f,np.float32,Nsubs)
                catalog['SubVmax'][istartSub:(istartSub+Nsubs)] = np.fromfile(f,np.float32,Nsubs)
                thisxyz = np.fromfile(f,np.float32,3*Nsubs)
                catalog['SubSpin'][istartSub:(istartSub+Nsubs),:] = np.reshape(thisxyz,[Nsubs,3])
                catalog['SubMostBoundID'][istartSub:(istartSub+Nsubs)] = np.fromfile(f,np.int64,Nsubs)
                catalog['SubHalfMass'][istartSub:(istartSub+Nsubs)] = np.fromfile(f,np.float32,Nsubs)

                istartSub += Nsubs
                istartGroup += Ngroups
                istartIDs += Nids
            
            f.close()

    return catalog


def getsubids(runid,snapnum,datadir=None,datascope=False,verbose=False):
    """Reads the subhalo particle IDs of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 128 to 511
        or as a length 3 tuple giving the 3-digit ID as (`X`,`Y`,`Z`)
        where `X` goes from 2 to 7 and `Y` and `Z` each go from 0 to 7.
        Ignored if `datadir` is set.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation `X_Y_Z`.
        Default is to read from the output of ``get_loc(runid,snapnum)``.
    datascope : boolean, optional
        If True, read from ``/datascope_path/indraX/X_Y_Z/`` (default False).
        Ignored if `datadir` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    subIDs : ndarray
        The IDs of each subhalo's constituent particles.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid,snapnum)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    idsfile = '{0}/postproc_{1:03d}/sub_ids_{1:03d}.'.format(datadir,snapnum)
    
    # open first file to get NTask just in case (but it should be 256) and TotNgroups
    f = open(idsfile+str(0),'rb')
    Ngroups, Nids, TotNgroups, NTask = np.fromfile(f, np.int32, 4)
    f.close()
    TotNsubs = getsubheader(runid,snapnum,datadir)

    if TotNsubs == 0: 
        if (run.X,run.Y,run.Z,snapnum) in MissingFiles_Problem:
            print("Missing data in {}_{}_{} snapshot {}: returning None".format(run.X,run.Y,run.Z,snapnum))
        else:
            print("No subhalos in {}_{}_{} snapshot {}: returning None".format(run.X,run.Y,run.Z,snapnum))

        return None

    else:
        TotSubids = 0
        # Get total number of IDs (including unbound particle IDs)
        for i in np.arange(0,NTask):
            f = open(idsfile+str(i),'rb')
            Ngroups, Nids, TotNgroups, NTask = np.fromfile(f, np.int32, 4)
            TotSubids += Nids
            f.close()
        subids = np.zeros(TotSubids,dtype=np.int64)
        
        # loop through numfiles (NTask)
        istart = 0
        for i in np.arange(0,NTask):
            f = open(idsfile+str(i),'rb')

            Ngroups, Nids, TotNgroups, NTask = np.fromfile(f, np.int32, 4)
            if Nids > 0:
                locIDs = np.fromfile(f,np.int64,Nids)
                # bitshift to remove the hash table info from the IDs
                subids[istart:(istart+Nids)] = np.bitwise_and(locIDs[:], (np.int64(1)<<34) - 1)
                istart = istart + Nids
            f.close()

    return subids-1


"""
FFT functions: 
fft_re,fft_im,kx,ky,kz are all arrays of the same shape: (129,129,65)
"""

def getfft(runid,tnum,datadir=None,datascope=False,verbose=False):
    """Reads the FFT output at one timestep of an Indra simulation.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 128 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where `X` goes from 2 to 7 and `Y` and `Z` each go from 0 to 7.
    tnum : int
        Which time step to read (0 to 504).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope_path/indraX/X_Y_Z/ (default False).
        Ignored if ``datadir`` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    fft_re, fft_im, time : tuple of (ndarray, ndarray, float)
        Real and Imaginary components of the Fourier modes of the 
        density field as arrays of shape [129,129,65], and the
        scalefactor of this time step.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid) #snapnum not specified, all files are on FileDB

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    filename = '{0}/FFT_DATA/FFT_128_{1:03d}.dat'.format(datadir,tnum)
    
    L = 128
    L2 = L//2

    # Check that files exist:
    try:
        with open(filename,'rb') as f:
            time = np.fromfile(f,np.float64,1)
            nsize = np.fromfile(f,np.int32,1)
            fft_re = np.fromfile(f,np.float32,nsize[0])
            fft_im = np.fromfile(f,np.float32,nsize[0])
            f.close()

            fft_re = np.reshape(fft_re,[L+1,L+1,L2+1])
            fft_im = np.reshape(fft_im,[L+1,L+1,L2+1])
            #print('a = {}'.format(time[0]))

            return fft_re,fft_im,time[0]
    except FileNotFoundError:
        # doesn't exist
        print('No file {}: returning None'.format(filename))
        return None,None,None


def getkvals(L=128):
    """Computes the k values that correspond the Fourier mode outputs.
    
    Parameters
    ----------
    L : int, optional
        The 1-dimensional grid size (defaults to 128).
    
    Returns
    -------
    kx, ky, kz : tuple of ndarrays
        The components of the *k*-values in the "upper-half sphere" of
        *k*-space (`kz` >= 0). Each array has size ``[L+1,L+1,L/2+1]``, the 
        same as the FFT outputs when `L`=128.
    """
    # define k's that correspond to fourier modes: (2*np.pi/boxsize)*np.array(x,y,z)
    # x = [-L/2:L/2], y = [-L/2:L/2], z = [0:L/2]
    
    L2 = L//2
    boxsize = 1000.

    kx = np.atleast_3d(np.expand_dims(np.arange(-L2,L2+1),axis=1))
    ky = np.atleast_3d(np.expand_dims(np.arange(-L2,L2+1),axis=0))
    kz = np.expand_dims(np.expand_dims(np.arange(0,L2+1),axis=0),axis=0)
    kx = np.broadcast_to(kx,(L+1,L+1,L2+1))*np.pi*2/boxsize
    ky = np.broadcast_to(ky,(L+1,L+1,L2+1))*np.pi*2/boxsize
    kz = np.broadcast_to(kz,(L+1,L+1,L2+1))*np.pi*2/boxsize
    return kx,ky,kz
