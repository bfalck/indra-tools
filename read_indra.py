"""
Indra reading functions using NumPy, tested in both Python 2 and Python 3

Inputs: 
- X, Y, and Z specify the Indra run, and are ignored if a datadir is specified
- datadir defaults to /datascope/indraX/X_Y_Z/
- snapnum goes from 0 to 63
- tnum goes from 0 to 504 for the FFT data

Usage:

--- Snapshots ---

header = getheader(X,Y,Z,snapnum,datadir=datadir)
pos = getpos(X,Y,Z,snapnum,datadir=datadir)
pos, vel, ids = getparticles(X,Y,Z,snapnum,datadir=datadir,sort=False)

--- Halo and subhalo data ---

TotNgroups, NTask = getfofheader(X,Y,Z,snapnum,datadir=datadir)
TotNsubs, NTask = getsubheader(X,Y,Z,snapnum,datadir=datadir)
groupLen, groupOffset = getfof(X,Y,Z,snapnum,datadir=None)
groupLen, groupOffset, groupids = getfofids(X,Y,Z,snapnum,datadir=None)
catalog = getsubcat(X,Y,Z,snapnum,datadir=datadir) # dictionary contains subLen and subOffset
subids = getsubids(X,Y,Z,snapnum,datadir=datadir)

--- FFT data ---

fft_re, fft_im, a = getfft(X,Y,Z,tnum,datadir=datadir)
kx, ky, kz = getkvals(L=128) # These are not read from file but are built and have the same shapes as fft_re and fft_im



Written by Bridget Falck, 2018

"""


import numpy as N


def readheader(f):

    if N.fromfile(f,N.int32,1) != 256:
        return None
    header={}
        
    npall = N.fromfile(f,N.int32,6) # DM only simulations - other np's and masses are 0
    header['np_file'] = npall[1] # number of particles in this file (needed to read positions)
    massall = N.fromfile(f,N.float64,6)
    header['mass'] = massall[1]
    header['time'],header['redshift'] = N.fromfile(f,N.float64,2)
    header['flag_sfr'],header['flag_feedback'] = N.fromfile(f,N.int32,2)
    nptotal = N.fromfile(f,N.int32,6)
    header['npart'] = nptotal[1]
    header['flag_cooling'],header['num_files'] = N.fromfile(f,N.int32,2)
    header['BoxSize'],header['omega0'],header['omegaL'],header['hubble'] = N.fromfile(f,N.float64,4)
    header['flag_stellarage'],header['flag_metals'],header['hashtabsize'] = N.fromfile(f,N.int32,3)
    # read rest of header_size + 2 filler integers:
    empty = N.fromfile(f,N.int32,23)
    return header

def getheader(X,Y,Z,snapnum,datadir=None):

    if (datadir == None): datadir = '/datascope/indra%s/%s_%s_%s/'%(str(X),str(X),str(Y),str(Z))

    sn = "%03d" % snapnum
    filename = datadir+'/snapdir_'+sn+'/snapshot_'+sn+'.'
    f = open(filename+str(0),'rb')

    header = readheader(f)
    f.close()
    
    return header

def readpos(f,npfile):
    
    thispos = N.fromfile(f,N.float32,3*npfile)
    thispos = N.reshape(thispos, [npfile, 3])
    
    return thispos

def readsnap(f,npfile):
    
    thispos = readpos(f,npfile)
    empty = N.fromfile(f,N.int32,2)
    thisvel = N.fromfile(f,N.float32,3*npfile)
    thisvel = N.reshape(thisvel, [npfile, 3])
    empty = N.fromfile(f,N.int32,2)
    thisID = N.fromfile(f,N.int64,npfile)
    
    return thispos,thisvel,thisID

def getpos(X,Y,Z,snapnum,datadir=None):
    
    if (datadir == None): datadir = '/datascope/indra%s/%s_%s_%s/'%(str(X),str(X),str(Y),str(Z))
    
    sn = "%03d" % snapnum
    filename = datadir+'/snapdir_'+sn+'/snapshot_'+sn+'.'

    NTask = 256
    nparticles = 1024
    
    # loop through files
    pos = N.empty((nparticles**3,3),N.float32)
    istart = 0
    for i in N.arange(0,NTask):
        f = open(filename+str(i),'rb')
        
        header = readheader(f)
        thispos = readpos(f,header['np_file'])
        
        f.close()
    
        pos[istart:(istart+header['np_file']),:] = thispos
        istart = istart + header['np_file']

    return pos

def getparticles(X,Y,Z,snapnum,datadir=None,sort=False):
    
    if (datadir == None): datadir = '/datascope/indra%s/%s_%s_%s/'%(str(X),str(X),str(Y),str(Z))

    sn = "%03d" % snapnum
    filename = datadir+'/snapdir_'+sn+'/snapshot_'+sn+'.'

    NTask = 256
    nparticles = 1024

    # loop through files
    pos = N.empty((nparticles**3,3),N.float32)
    vel = N.empty((nparticles**3,3),N.float32)
    ids = N.empty((nparticles**3),N.int64)
    istart = 0
    for i in N.arange(0,NTask):
        f = open(filename+str(i),'rb')
        header = readheader(f)
        thispos,thisvel,thisID = readsnap(f,header['np_file'])
        f.close()

        thisID -= 1 # make IDs start at 0 to use to index pos and vel
        
        if sort:
            pos[thisID] = thispos
            vel[thisID] = thisvel
        else:
            pos[istart:(istart+header['np_file']),:] = thispos
            vel[istart:(istart+header['np_file']),:] = thisvel
            ids[istart:(istart+header['np_file'])] = thisID
    
        istart = istart + header['np_file']
        
    if sort: ids = N.arange(header['npart'])

    return pos,vel,ids


"""
Halo and subhalo functions
"""

def getfofheader(X,Y,Z,snapnum,datadir=None):
    
    if (datadir == None): datadir = '/datascope/indra%s/%s_%s_%s/'%(str(X),str(X),str(Y),str(Z))

    sn = "%03d" % snapnum
    tabfile = datadir+'snapdir_'+sn+'/group_tab_'+sn+'.'

    f = open(tabfile+str(0),'rb')
    Ngroups, Nids, TotNgroups, NTask = N.fromfile(f, N.int32, 4)
    f.close()
    
    return TotNgroups,NTask

def getfof(X,Y,Z,snapnum,datadir=None):
    
    if (datadir == None): datadir = '/datascope/indra%s/%s_%s_%s/'%(str(X),str(X),str(Y),str(Z))

    sn = "%03d" % snapnum
    tabfile = datadir+'snapdir_'+sn+'/group_tab_'+sn+'.'

    # loop through NTask files (could read NTask from header...)
#    NTask = 256
    TotNgroups,NTask = getfofheader(X,Y,Z,snapnum,datadir)
    # Don't read if no groups...
    if TotNgroups == 0: return 0,0
    else:
        groupLen = N.zeros(TotNgroups,dtype=N.int32)
        groupOffset = N.zeros(TotNgroups,dtype=N.int32)

        istartGroup = 0
        istartIDs = 0
        for i in N.arange(0,NTask):
            f = open(tabfile+str(i), 'rb')
    
            Ngroups, Nids, TotNgroups, NTask = N.fromfile(f, N.int32, 4)
            if Ngroups > 0:
                locLen = N.fromfile(f,N.int32,Ngroups)
                locOffset = N.fromfile(f,N.int32,Ngroups)
                groupLen[istartGroup:(istartGroup+Ngroups)] = locLen
                groupOffset[istartGroup:(istartGroup+Ngroups)] = locOffset+istartIDs
                istartGroup += Ngroups
                istartIDs += Nids
#            else: print('No groups in file %d' % i)
            f.close()
    
    return groupLen,groupOffset

def getfofids(X,Y,Z,snapnum,datadir=None):

    if (datadir == None): datadir = '/datascope/indra%s/%s_%s_%s/'%(str(X),str(X),str(Y),str(Z))

    sn = "%03d" % snapnum
    idsfile = datadir+'snapdir_'+sn+'/group_ids_'+sn+'.'
 
    # loop through NTask files
#    NTask = 256
    TotNgroups,NTask = getfofheader(X,Y,Z,snapnum,datadir)
    # Don't read if no groups...
    if TotNgroups == 0: return 0,0
    else:
        groupLen,groupOffset = getfof(X,Y,Z,snapnum,datadir)
        TotNids = N.sum(groupLen,dtype=N.int64)
        groupids = N.zeros(TotNids,dtype=N.int64)

        istart = 0
        for i in N.arange(0,NTask):
            f = open(idsfile+str(i),'rb')

            Ngroups, Nids, TotNgroups, NTask = N.fromfile(f, N.int32, 4)
            if Nids > 0:
                locIDs = N.fromfile(f,N.int64,Nids)
                # bitshift to remove the hash table info from the IDs
                groupids[istart:(istart+Nids)] = N.bitwise_and(locIDs[:], (N.int64(1)<<34) - 1)
                istart += Nids
#            else: print('No groups in file %d' % i)
            f.close()
    
    return groupLen,groupOffset,groupids-1

def getsubheader(X,Y,Z,snapnum,datadir=None):

    if (datadir == None): datadir = '/datascope/indra%s/%s_%s_%s/'%(str(X),str(X),str(Y),str(Z))

    sn = "%03d" % snapnum
    tabfile = datadir+'postproc_'+sn+'/sub_tab_'+sn+'.'

    TotNsubs = 0
    f = open(tabfile+str(0),'rb')
    Ngroups,Nids,TotNgroups,NTask = N.fromfile(f,N.int32,4)
    f.close()
    for i in N.arange(0,NTask):
        f = open(tabfile+str(i),'rb')
#        print('opening file '+tabfile+str(i))
        Ngroups, Nids, TotNgroups, NTask, Nsubs = N.fromfile(f, N.int32, 5)
        TotNsubs += Nsubs
        f.close()

    return TotNsubs,NTask

def getsubcat(X,Y,Z,snapnum,datadir=None):
    
    if (datadir == None): datadir = '/datascope/indra%s/%s_%s_%s/'%(str(X),str(X),str(Y),str(Z))

    sn = "%03d" % snapnum
    tabfile = datadir+'postproc_'+sn+'/sub_tab_'+sn+'.'
    
    TotNgroups,NTask= getfofheader(X,Y,Z,snapnum,datadir)
    TotNsubs,NTask = getsubheader(X,Y,Z,snapnum,datadir)
    if TotNsubs == 0: return None
    else:
        # Initialize data containers. Group into structured array?? Dict??
        catalog = {}
        catalog['NsubPerHalo'] = N.zeros(TotNgroups,dtype=N.int32)
        catalog['FirstSubOfHalo'] = N.zeros(TotNgroups,dtype=N.int32) # file specific!
    
        catalog['subLen'] = N.zeros(TotNsubs,dtype=N.int32)
        catalog['subOffset'] = N.zeros(TotNsubs,dtype=N.int32) # file specific!
        catalog['subParentHalo'] = N.zeros(TotNsubs,dtype=N.int32) # file specific!
    
        catalog['M200mean'] = N.zeros(TotNgroups,dtype=N.float32)
        catalog['R200mean'] = N.zeros(TotNgroups,dtype=N.float32)
        catalog['M200crit'] = N.zeros(TotNgroups,dtype=N.float32)
        catalog['R200crit'] = N.zeros(TotNgroups,dtype=N.float32)
        catalog['M200tophat'] = N.zeros(TotNgroups,dtype=N.float32)
        catalog['R200tophat'] = N.zeros(TotNgroups,dtype=N.float32)
    
        catalog['SubPos'] = N.zeros((TotNsubs,3),dtype=N.float32)
        catalog['SubVel'] = N.zeros((TotNsubs,3),dtype=N.float32)
        catalog['SubVelDisp'] = N.zeros(TotNsubs,dtype=N.float32)
        catalog['SubVmax'] = N.zeros(TotNsubs,dtype=N.float32)
        catalog['SubSpin'] = N.zeros((TotNsubs,3),dtype=N.float32)
        catalog['SubMostBoundID'] = N.zeros((TotNsubs,2),dtype=N.int32)
        catalog['SubHalfMass'] = N.zeros(TotNsubs,dtype=N.float32)
        
        # loop over numfiles (NTask)
        istartSub = 0
        istartGroup = 0
        istartIDs = 0
        for i in N.arange(0,NTask):
            f = open(tabfile+str(i),'rb')
            Ngroups, Nids, TotNgroups, NTask, Nsubs = N.fromfile(f, N.int32, 5)
            
            if Nsubs > 0:
            # Read catalog: indexes need to include offset from previous files (istarts)
                catalog['NsubPerHalo'][istartGroup:(istartGroup+Ngroups)] = N.fromfile(f,N.int32,Ngroups)
                catalog['FirstSubOfHalo'][istartGroup:(istartGroup+Ngroups)] = N.fromfile(f,N.int32,Ngroups)+istartSub
       
                catalog['subLen'][istartSub:(istartSub+Nsubs)] = N.fromfile(f,N.int32,Nsubs)
                catalog['subOffset'][istartSub:(istartSub+Nsubs)] = N.fromfile(f,N.int32,Nsubs)+istartIDs
                catalog['subParentHalo'][istartSub:(istartSub+Nsubs)] = N.fromfile(f,N.int32,Nsubs)+istartGroup
                
                catalog['M200mean'][istartGroup:(istartGroup+Ngroups)] = N.fromfile(f,N.float32,Ngroups)
                catalog['R200mean'][istartGroup:(istartGroup+Ngroups)] = N.fromfile(f,N.float32,Ngroups)
                catalog['M200crit'][istartGroup:(istartGroup+Ngroups)] = N.fromfile(f,N.float32,Ngroups)
                catalog['R200crit'][istartGroup:(istartGroup+Ngroups)] = N.fromfile(f,N.float32,Ngroups)
                catalog['M200tophat'][istartGroup:(istartGroup+Ngroups)] = N.fromfile(f,N.float32,Ngroups)
                catalog['R200tophat'][istartGroup:(istartGroup+Ngroups)] = N.fromfile(f,N.float32,Ngroups)

                thisxyz = N.fromfile(f,N.float32,3*Nsubs)
                catalog['SubPos'][istartSub:(istartSub+Nsubs),:] = N.reshape(thisxyz,[Nsubs,3])
                thisxyz = N.fromfile(f,N.float32,3*Nsubs)
                catalog['SubVel'][istartSub:(istartSub+Nsubs),:] = N.reshape(thisxyz,[Nsubs,3])
                catalog['SubVelDisp'][istartSub:(istartSub+Nsubs)] = N.fromfile(f,N.float32,Nsubs)
                catalog['SubVmax'][istartSub:(istartSub+Nsubs)] = N.fromfile(f,N.float32,Nsubs)
                thisxyz = N.fromfile(f,N.float32,3*Nsubs)
                catalog['SubSpin'][istartSub:(istartSub+Nsubs),:] = N.reshape(thisxyz,[Nsubs,3])
                catalog['SubMostBoundID'][istartSub:(istartSub+Nsubs),:] = N.reshape(N.fromfile(f,N.int32,2*Nsubs),[Nsubs,2])
                catalog['SubHalfMass'][istartSub:(istartSub+Nsubs)] = N.fromfile(f,N.float32,Nsubs)

                istartSub += Nsubs
                istartGroup += Ngroups
                istartIDs += Nids
            
            f.close()

    return catalog


def getsubids(X,Y,Z,snapnum,datadir=None):
    
    if (datadir == None): datadir = '/datascope/indra%s/%s_%s_%s/'%(str(X),str(X),str(Y),str(Z))

    sn = "%03d" % snapnum
    idsfile = datadir+'postproc_'+sn+'/sub_ids_'+sn+'.'
    
    TotNsubs,NTask = getsubheader(X,Y,Z,snapnum,datadir)
    if TotNsubs == 0: return None
    else:
        TotSubids = 0
        # Get total number of IDs (including unbound particle IDs)
        for i in N.arange(0,NTask):
            f = open(idsfile+str(i),'rb')
            Ngroups, Nids, TotNgroups, NTask = N.fromfile(f, N.int32, 4)
            TotSubids += Nids
            f.close()
        subids = N.zeros(TotSubids,dtype=N.int64)
        
        # loop through numfiles (NTask)
        istart = 0
        for i in N.arange(0,NTask):
            f = open(idsfile+str(i),'rb')

            Ngroups, Nids, TotNgroups, NTask = N.fromfile(f, N.int32, 4)
            if Nids > 0:
                locIDs = N.fromfile(f,N.int64,Nids)
                # bitshift to remove the hash table info from the IDs
                subids[istart:(istart+Nids)] = N.bitwise_and(locIDs[:], (N.int64(1)<<34) - 1)
                istart = istart + Nids
            f.close()

    return subids-1


"""
FFT functions: 
fft_re,fft_im,kx,ky,kz are all arrays of the same shape: (129,129,65)
"""

def getfft(X,Y,Z,tnum,datadir=None):
    # Read FFT data in given time slice (0 to 504)
    
    if (datadir == None): datadir = '/datascope/indra%s/%s_%s_%s/'%(str(X),str(X),str(Y),str(Z))
    
    tstr = "%03d" % tnum
    filename = datadir+'FFT_DATA/FFT_128_%s.dat' % tstr  

    L = 128
    L2 = L//2

    f = open(filename,'rb')
    time = N.fromfile(f,N.float64,1)
    nsize = N.fromfile(f,N.int32,1)
    fft_re = N.fromfile(f,N.float32,nsize[0])
    fft_im = N.fromfile(f,N.float32,nsize[0])
    f.close()

    fft_re = N.reshape(fft_re,[L+1,L+1,L2+1])
    fft_im = N.reshape(fft_im,[L+1,L+1,L2+1])
    #print 'a = %f' % time[0]

    return fft_re,fft_im,time[0]

def getkvals(L=None):
    # define k's that correspond to fourier modes: (2*N.pi/boxsize)*N.array(x,y,z)
    # x = [-L/2:L/2], y = [-L/2,L/2], z = [0,L/2]
    # L defaults to 128, as in FFT_DATA output
    
    if (L == None): L = 128
    L2 = L//2
    boxsize = 1000.

    kx = N.atleast_3d(N.expand_dims(N.arange(-L2,L2+1),axis=1))
    ky = N.atleast_3d(N.expand_dims(N.arange(-L2,L2+1),axis=0))
    kz = N.expand_dims(N.expand_dims(N.arange(0,L2+1),axis=0),axis=0)
    kx = N.broadcast_to(kx,(L+1,L+1,L2+1))*N.pi*2/boxsize
    ky = N.broadcast_to(ky,(L+1,L+1,L2+1))*N.pi*2/boxsize
    kz = N.broadcast_to(kz,(L+1,L+1,L2+1))*N.pi*2/boxsize
    return kx,ky,kz