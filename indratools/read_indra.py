"""
Reading functions for the Indra suite of simulations hosted on the SciServer.

Written by Bridget Falck, 2018-2020


TODO: add usage examples
TODO: add/update error handling


Inputs: 
- ``runid`` specifies the Indra run, and is ignored if a ``datadir`` is specified; 
    it is either an integer from 0 to 511 or a tuple containing (``X``,``Y``,``Z``)
- ``datadir`` defaults to the FileDB location of run ``X_Y_Z``. If ``datadir`` is not 
    set and ``datascope=True``, the ``datadir`` will be the datascope location of 
    run ``X_Y_Z``, e.g. ``/datascope/indraX/X_Y_Z/``.
- ``snapnum`` goes from 0 to 63
- ``tnum`` goes from 0 to 504 for the FFT data


Methods
-------

--- Snapshots ---

getheader(runid,snapnum,datadir=None,datascope=False,verbose=False)
    Reads and returns a dictionary of header parameters.
getpos(runid,snapnum,datadir=None,datascope=False,verbose=False)
    Reads and returns an array of particle positions of shape [1024**3,3]
    from one full snapshot.
getparticles(runid,snapnum,datadir=None,datascope=False,sort=False,verbose=False)
    Reads particle positions, velocities, and IDs: pos, vel, ids = getparticles(...).
    pos and vel are arrays of shape [1024**3,3], and ids have shape [1024**3].

--- Halo and subhalo data ---

getfofheader(runid,snapnum,datadir=None,datascope=False,verbose=False)
    Reads the header of the FOF data and returns the total number of FOF 
    groups contained in all files at this snapshot.
getsubheader(runid,snapnum,datadir=None,datascope=False,verbose=False)
    Reads all headers of the SUBFIND files and returns the total number
    of subhalos in all files at this snapshot.
getfof(runid,snapnum,datadir=None,datascope=False,verbose=False)
    Reads the number of particles in each FOF group and the Offset array needed to
    index the IDs of the member particles: groupLen, groupOffset = getfof(...)
getfofids(runid,snapnum,datadir=None,datascope=False,verbose=False)
    Reads the groupLen, groupOffset, and particle ID arrays for the FOF groups.
getsubcat(runid,snapnum,datadir=None,datascope=False,verbose=False)
    Reads the SUBFIND subhalo catalogs and returns a dictionary of values. Some halo 
    properties are defined for each FOF group, and some for each subhalo. The
    dictionary contains the subLen and subOffset arrays needed to index the IDs
    of the subhalo member particles (as for the FOF groups).
getsubids(runid,snapnum,datadir=None,datascope=False,verbose=False)
    Reads and returns the IDs of the particles in each subhalo.

--- FFT data ---

getfft(runid,tnum,datadir=None,datascope=False,verbose=False)
    Reads and returns the real and imaginary components of the Fourier-space density
    grid at output ``tnum`` and returns the scalefactor of this ``tnum``.
getkvals(L=128) 
    Builds and returns the x, y, and z components of the k-vectors associated
    with the FFT data. Each are arrays with the same shapes as ``fft_re`` and ``fft_im``.
"""


import numpy as np


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
            self.X, self.Y, self.Z = np.unravel_index(runid,(8,8,8))
        elif isinstance(runid, tuple):
            self.X, self.Y, self.Z = runid
            self.num = np.ravel_multi_index(runid,(8,8,8))


def get_loc(runid):
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


def readheader(f):
    """Utility function for reading Indra snapshot files"""

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


def getheader(runid,snapnum,datadir=None,datascope=False,verbose=False):
    """Reads the Gadget header of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope/indraX/X_Y_Z/ (default False).
        Overwritten if ``datascope`` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    header : dict
        Dictionary object containing header information.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = '/datascope/indra{0}/{0}_{1}_{2}/'.format(run.X,run.Y,run.Z)
        else:
            datadir = get_loc(runid)

    if (verbose == True): print('Reading from {}'.format(datadir))

    filename = '{0}/snapdir_{1:03d}/snapshot_{1:03d}.'.format(datadir,snapnum)
    f = open(filename+str(0),'rb')

    header = readheader(f)
    f.close()
    
    return header

def readpos(f,npfile):
    """Utility function for ``getpos()``"""
    
    thispos = np.fromfile(f,np.float32,3*npfile)
    thispos = np.reshape(thispos, [npfile, 3])
    
    return thispos

def readsnap(f,npfile):
    """Utility function for ``getparticles()``"""
    
    thispos = readpos(f,npfile)
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
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope/indraX/X_Y_Z/ (default False).
        Overwritten if ``datascope`` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    pos : ndarray
        Particle positions (in Mpc/h) in array of shape [1024**3,3]
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = '/datascope/indra{0}/{0}_{1}_{2}/'.format(run.X,run.Y,run.Z)
        else:
            datadir = get_loc(runid)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    filename = '{0}/snapdir_{1:03d}/snapshot_{1:03d}.'.format(datadir,snapnum)

    NTask = 256
    nparticles = 1024
    
    # loop through files
    pos = np.empty((nparticles**3,3),np.float32)
    istart = 0
    for i in np.arange(0,NTask):
        f = open(filename+str(i),'rb')
        
        header = readheader(f)
        thispos = readpos(f,header['np_file'])
        
        f.close()
    
        pos[istart:(istart+header['np_file']),:] = thispos
        istart = istart + header['np_file']

    return pos

def getparticles(runid,snapnum,datadir=None,datascope=False,sort=False,verbose=False):
    """Reads the particle positions, velocities, and IDs of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope/indraX/X_Y_Z/ (default False).
        Overwritten if ``datascope`` is set.
    sort : boolean, optional
        If True, sort the positions and velocities by the particle IDs.
        The sorting is done on-the-fly to minimize reading time.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    pos, vel, ids : tuple of ndarrays
        Particle positions (in Mpc/h) and velocities (in km/s) in arrays of 
        shape [1024**3,3] and IDs in array of shape [1024**3]. If ``sort==True``,
        ID array is equivalent to ``np.arange(1024**3)``.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = '/datascope/indra{0}/{0}_{1}_{2}/'.format(run.X,run.Y,run.Z)
        else:
            datadir = get_loc(runid)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    filename = '{0}/snapdir_{1:03d}/snapshot_{1:03d}.'.format(datadir,snapnum)

    NTask = 256
    nparticles = 1024

    # loop through files
    pos = np.empty((nparticles**3,3),np.float32)
    vel = np.empty((nparticles**3,3),np.float32)
    ids = np.empty((nparticles**3),np.int64)
    istart = 0
    for i in np.arange(0,NTask):
        f = open(filename+str(i),'rb')
        header = readheader(f)
        thispos,thisvel,thisID = readsnap(f,header['np_file'])
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
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope/indraX/X_Y_Z/ (default False).
        Overwritten if ``datascope`` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    TotNgroups : int
        Total number of FOF groups in this snapshot over all 256 files.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = '/datascope/indra{0}/{0}_{1}_{2}/'.format(run.X,run.Y,run.Z)
        else:
            datadir = get_loc(runid)

    if (verbose == True): print('Reading from {}'.format(datadir))

    tabfile = '{0}/snapdir_{1:03d}/group_tab_{1:03d}.'.format(datadir,snapnum)

    f = open(tabfile+str(0),'rb')
    Ngroups, Nids, TotNgroups, NTask = np.fromfile(f, np.int32, 4)
    f.close()
    
    return TotNgroups

def getfof(runid,snapnum,datadir=None,datascope=False,verbose=False):
    """Reads the FOF group info of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope/indraX/X_Y_Z/ (default False).
        Overwritten if ``datascope`` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    groupLen, groupOffset : tuple of ndarrays
        Number of particles in each FOF group, and Offset array needed
        in order to index the IDs of the group's constituent particles.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = '/datascope/indra{0}/{0}_{1}_{2}/'.format(run.X,run.Y,run.Z)
        else:
            datadir = get_loc(runid)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    tabfile = '{0}/snapdir_{1:03d}/group_tab_{1:03d}.'.format(datadir,snapnum)

    # loop through NTask files
    NTask = 256
    TotNgroups = getfofheader(runid,snapnum,datadir)
    # Don't read if no groups...
    if TotNgroups == 0: return 0,0
    else:
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
#            else: print('No groups in file %d' % i)
            f.close()
    
    return groupLen,groupOffset

def getfofids(runid,snapnum,datadir=None,datascope=False,verbose=False):
    """Reads the FOF group info (including particle IDs) of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope/indraX/X_Y_Z/ (default False).
        Overwritten if ``datascope`` is set.
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
        if (datascope == True): datadir = '/datascope/indra{0}/{0}_{1}_{2}/'.format(run.X,run.Y,run.Z)
        else:
            datadir = get_loc(runid)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    idsfile = '{0}/snapdir_{1:03d}/group_ids_{1:03d}.'.format(datadir,snapnum)
 
    # loop through NTask files
    NTask = 256
    TotNgroups = getfofheader(runid,snapnum,datadir)
    # Don't read if no groups...
    if TotNgroups == 0: return 0,0
    else:
        groupLen,groupOffset = getfof(runid,snapnum,datadir)
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

def getsubheader(runid,snapnum,datadir=None,datascope=False,verbose=False):
    """Reads the SUBFIND header (total number of subhalos) of an Indra snapshot.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope/indraX/X_Y_Z/ (default False).
        Overwritten if ``datascope`` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    TotNsubs : int
        Total number of subhalos in this snapshot over all 256 files.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = '/datascope/indra{0}/{0}_{1}_{2}/'.format(run.X,run.Y,run.Z)
        else:
            datadir = get_loc(runid)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    tabfile = '{0}/postproc_{1:03d}/sub_tab_{1:03d}.'.format(datadir,snapnum)

    TotNsubs = 0
    f = open(tabfile+str(0),'rb')
    Ngroups,Nids,TotNgroups,NTask = np.fromfile(f,np.int32,4)
    f.close()
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
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope/indraX/X_Y_Z/ (default False).
        Overwritten if ``datascope`` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    catalog : dictionary
        Dictionary of halo and subhalo properties and indexing information.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = '/datascope/indra{0}/{0}_{1}_{2}/'.format(run.X,run.Y,run.Z)
        else:
            datadir = get_loc(runid)

#    print('Reading from {}'.format(datadir))
    
    tabfile = '{0}/postproc_{1:03d}/sub_tab_{1:03d}.'.format(datadir,snapnum)
    
    NTask = 256
    TotNgroups= getfofheader(runid,snapnum,datadir)
    TotNsubs = getsubheader(runid,snapnum,datadir)
    if TotNsubs == 0: return None
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
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int
        Which snapshot to read (0 to 63).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope/indraX/X_Y_Z/ (default False).
        Overwritten if ``datascope`` is set.
    verbose : boolean, optional
        If True, print reading statements (default False).
    
    Returns
    -------
    subIDs : ndarray
        The IDs of each subhalo's constituent particles.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = '/datascope/indra{0}/{0}_{1}_{2}/'.format(run.X,run.Y,run.Z)
        else:
            datadir = get_loc(runid)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    idsfile = '{0}/postproc_{1:03d}/sub_ids_{1:03d}.'.format(datadir,snapnum)
    
    NTask = 256
    TotNsubs = getsubheader(runid,snapnum,datadir)
    if TotNsubs == 0: return None
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
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    tnum : int
        Which time step to read (0 to 504).
    datadir : string, optional
        If set, specify full path of directory containing simulation X_Y_Z.
        Default is to read from the output of ``get_loc(runid)``.
    datascope : boolean, optional
        If True, read from /datascope/indraX/X_Y_Z/ (default False).
        Overwritten if ``datascope`` is set.
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
        if (datascope == True): datadir = '/datascope/indra{0}/{0}_{1}_{2}/'.format(run.X,run.Y,run.Z)
        else:
            datadir = get_loc(runid)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    tabfile = '{0}/FFT_DATA/FFT_128_{1:03d}.dat'.format(datadir,tnum)
    
    L = 128
    L2 = L//2

    f = open(filename,'rb')
    time = np.fromfile(f,np.float64,1)
    nsize = np.fromfile(f,np.int32,1)
    fft_re = np.fromfile(f,np.float32,nsize[0])
    fft_im = np.fromfile(f,np.float32,nsize[0])
    f.close()

    fft_re = np.reshape(fft_re,[L+1,L+1,L2+1])
    fft_im = np.reshape(fft_im,[L+1,L+1,L2+1])
    #print('a = {}'.format(time[0]))

    return fft_re,fft_im,time[0]

def getkvals(L=128):
    """Computes the k values that correspond the Fourier mode outputs.
    
    Parameters
    ----------
    L : int, optional
        The 1-dimensional grid size (defaults to 128).
    
    Returns
    -------
    kx, ky, kz : tuple of ndarrays
        The components of the k-values in the "upper-half sphere" of
        k-space (kz >= 0). Each array has size [L+1,L+1,L/2+1], the 
        same as the FFT outputs when ``L=128``.
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
