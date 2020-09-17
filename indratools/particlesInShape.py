"""
Functions to efficiently read particles in Indra simulations that
are within geometrical shapes defined by Shape3D objects.

Written by Bridget Falck and Gerard Lemson


TO DO: Decide whether to include option to return PH cell data (results of query)


Methods
-------

get_run_num(x,y,z)
    Helper function to return the raveled run number (0 to 511) from the unraveled
    x, y, z identifiers (each 0 to 7).

"""

from .utils import *
import numpy as np
import SciServer.CasJobs as cj
import pandas


ds_basedir = '/home/idies/workspace/indra_dss/'


def _readParticles(file,shift,gr,shape,getvel,getIDs):
    """Does the actual file reading according to file and PH cell info in gr,
    then calls shape.contained() to select particles inside the shape."""
    npart = sum(gr['ixcount'].values)
    nalloc = npart*3
    parr = np.zeros(nalloc,dtype=np.float32)
    if getvel: varr = np.zeros(nalloc,dtype=np.float32)
    if getIDs: idarr = np.zeros(npart,dtype=np.int64)
    ncount = 0; idcount = 0
    # open file
    with open(file,"rb") as f:
        # if velocities or IDs requested, need to read number of particles in this file
        if getvel or getIDs:
            if np.fromfile(f,np.int32,1) != 256:
                return None
            npall = np.fromfile(f,np.int32,6)
            npfile = npall[1]
            # then calculate velshift = 3*npfile*sizeof(float32) PLUS fortran buffer: 2*int32
            velshift = 3*npfile*np.dtype(np.float32).itemsize + 2*np.dtype(np.int32).itemsize
            idshift = 2*velshift+4+256+4+4 # all pos, all vel, + header and empty bits
    # for gr tuples step through file
        for g in gr.itertuples():
            ix=g[4] # ibstart
            nb=g[5] # ibcount
            n=int(nb/4) # = g[3]*3 = 3*ixcount
            f.seek(ix) # for velocities or IDs, ix is larger by set amount
            chunk=f.read(nb)
            parr[ncount:(ncount+n)] = np.frombuffer(chunk,dtype=np.float32,count=n)
            if getvel:
                f.seek(ix+velshift)
                chunk=f.read(nb)
                varr[ncount:(ncount+n)] = np.frombuffer(chunk,dtype=np.float32,count=n)
            if getIDs:
                f.seek(8*g[2]+idshift) # ix skips 3*ixstart*float32, we need to skip ixstart*int64
                n_id = g[3]
                nb = 8*n_id
                chunk = f.read(nb)
                idarr[idcount:(idcount+n_id)] = np.frombuffer(chunk,dtype=np.int64,count=n_id)
                idcount += n_id
            ncount += n
    parr = parr.reshape(npart,3)
    inside = shape.contained(parr) # box defaults to 1000

    if getvel: varr = varr.reshape(npart,3)
    if getIDs: idarr -= 1 # want IDs to start at 0
    if getvel and getIDs:
        return parr[inside]+shift,varr[inside],idarr[inside]
    elif getvel:
        return parr[inside]+shift,varr[inside]
    elif getIDs:
        return parr[inside]+shift,idarr[inside]
    else:
        return parr[inside]+shift

def _retrieveParticles(df,shape,getvel,getIDs):
    """Loops through the file and PH cell info in df to call _readParticles,
    then re-organizes data into arrays not grouped by PH cell."""
    pos = {}
    if getvel: vel = {}
    if getIDs: ids = {}
    ntot = 0
    for k,gr in df.groupby(['loc','shiftx','shifty','shiftz']):
        particles = _readParticles(k[0],np.array(k[1:]),gr,shape,getvel,getIDs) 
        if particles is not None:
            if getvel and getIDs:
                pos[k] = particles[0]
                vel[k] = particles[1]
                ids[k] = particles[2]
            elif getvel:
                pos[k] = particles[0]
                vel[k] = particles[1]
            elif getIDs:
                pos[k] = particles[0]
                ids[k] = particles[1]
            else:
                pos[k] = particles
            if getvel or getIDs:
                ntot+=len(particles[0])
            else: ntot+=len(particles)
        else:
            print('No particles found in {0}'.format(k[0])) # maybe don't want to print a message... this occurs when PH cell in shape but individual particles in that cell aren't
    
    # re-organize particles (currently grouped by PH cell)
    part = {'NumParticles':ntot,
         'x':np.zeros(ntot),'y':np.zeros(ntot),'z':np.zeros(ntot)}
    if getvel:
        part['vx'] = np.zeros(ntot)
        part['vy'] = np.zeros(ntot)
        part['vz'] = np.zeros(ntot)
    if getIDs:
        part['ids'] = np.zeros(ntot,dtype=np.int64)
    count = 0
    # pos.keys() are same as vel.keys and ids.keys
    for k in pos.keys():
        npart = len(pos[k])
        part['x'][count:(count+npart)] = pos[k][:,0]
        part['y'][count:(count+npart)] = pos[k][:,1]
        part['z'][count:(count+npart)] = pos[k][:,2]
        if getvel:
            part['vx'][count:(count+npart)] = vel[k][:,0]
            part['vy'][count:(count+npart)] = vel[k][:,1]
            part['vz'][count:(count+npart)] = vel[k][:,2]
        if getIDs:
            part['ids'][count:(count+npart)] = ids[k]
        count+=npart

    return part

def particlesInShape(runid,snapnum,shape,getvel=False,getIDs=False,datadir=None,datascope=False,verbose=False):
    """
    Returns particles from given simulation and snapshot that
    are within shape. If shape crosses boundary, returned positions
    are wrapped around, so they can be negative.
    
    Parameters
    ----------
    runid : int or tuple
        Specifies the Indra run either as an integer from 0 to 511
        or as a length 3 tuple giving the 3-digit ID as (X,Y,Z)
        where X, Y, and Z each go from 0 to 7.
    snapnum : int
        Which snapshot to read (0 to 63).
    shape : Shape3D object
        A Box, Sphere, Cone, or ConeSegment object as defined in Shape3D.py
    getvel : bool (default False)
        Whether to read and return velocities in addition to positions.
    getIDs : bool (default False)
        Whether to read and return particle IDs in addition to positions.
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
    dict of {'NumParticles': int; 'x','y','z' [,'vx','vy','vz','ids'] : ndarrays}
        Requested information about the particles in the shape. Velocities and
        IDs only returned if getvel or getIDs are set to True. All arrays have
        length 'NumParticles'.
    """

    run = Run(runid)

    if (datadir == None): 
        if (datascope == True): datadir = f'/{ds_basedir}/indra{run.X}/{run.X}_{run.Y}_{run.Z}/'
        else:
            datadir = get_loc(runid)

    if (verbose == True): print('Reading from {}'.format(datadir))
    
    simId = f'{run.X}{run.Y}{run.Z}'
    
    sql=f"""
        select sciserverLocation as loc
        ,      ixstart as ixstart
        ,      ixcount as ixcount 
        ,      4+256+4+4+4*3*ixstart as ibstart
        ,      4*3*ixcount as ibcount 
        ,      shiftx
        ,      shifty
        ,      shiftz
          from fGetIndraBins({simId},{snapnum},'{str(shape)}') 
         order by sciserverLocation, shiftx, shifty, shiftz
    """
    _df = cj.executeQuery(sql,"Indra")
    part = _retrieveParticles(_df,shape,getvel,getIDs)

    return part
