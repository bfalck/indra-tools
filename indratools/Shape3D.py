"""
Object definitions for 3D shapes in N-body simulations. 

Written by Bridget Falck and Gerard Lemson. 
Based on Search3D library written by Gerard Lemson and Tamas Budavari.

Classes currently defined for:

- Sphere
- Box
- Cone
- ConeSegment

"""

import numpy as np

class Sphere():
    """
    Represents a spherical region. Can check for inclusion of a point and serialize itself in the 
    format accepted by the Spatial3D library's parser.
    
    Parameters
    ----------
    (x, y, z) : float
        The coordinates of the center of the sphere in simulation units.
    r : float
        The radius of the sphere in simulation units.
    
    """
    def __init__(self,x,y,z,r):
        
        if (r <= 0):
            raise ValueError("The radius of a sphere must be greater than 0.")
        self.x=x
        self.y=y
        self.z=z
        self.pos=[x,y,z]
        self.r=r
        self.r2=r*r
    
    def __str__(self):
        return f"SPHERE[{self.x:f},{self.y:f},{self.z:f},{self.r:f}]"
    
    def contained(self,pos,box=1000.):
        """
        Return a boolean index that will select particles in the shape.

        Parameters
        ----------
        pos : ndarray 
            Particle coordinate array of shape ``[number_of_particles,3]``
        box : float (default 1000)
            Size of simulation in Mpc/*h*.
            
        Returns
        -------
        boolean 
            Array of booleans giving whether each particle is within the shape.
        """
        hbox = box/2
        a = np.abs(pos-self.pos)
        a[a>hbox] -= box
        d = np.sum(a**2,axis=1)
        inside = (d <= self.r2)
        return inside


class Box():
    """
    Represents a 3-dimensional rectangular box parallel to the coordinate axes. Can check for 
    inclusion of a point and serialize itself in the format accepted by the Spatial3D library's parser.
    
    Parameters
    ----------
    (xmin, ymin, zmin) : float
        The coordinates of the lower-left corner of the box in simulation units.
        May be < 0 for a box that spans the periodic boundary.
    (xmax, ymax, zmax) : float
        The coordinates of the upper-right corner of the box in simulation units.
        May be greater than the simulation's length for a box that spans the periodic boundary.

    """
    def __init__(self,xmin,ymin,zmin,xmax,ymax,zmax):

        # Defined box must be in wrapped-around coordinates, i.e. corners can be < 0 or > box_size
        if ((xmin > xmax) or (ymin > ymax) or (zmin > zmax)):
            raise ValueError("Left-lower point of box must be left-lower from right-upper point; add or subtract L to span periodic boundary.")
        
        self.xmin = xmin
        self.ymin = ymin
        self.zmin = zmin
        self.xmax = xmax
        self.ymax = ymax
        self.zmax = zmax
        self.posmin = [xmin,ymin,zmin]
        self.posmax = [xmax,ymax,zmax]
        
    def __str__(self):
        return f"BOX[{self.xmin:f},{self.ymin:f},{self.zmin:f},{self.xmax:f},{self.ymax:f},{self.zmax:f}]"
    
    def contained(self,pos,box=1000):
        """
        Return a boolean index that will select particles in the shape.

        Parameters
        ----------
        pos : ndarray 
            Particle coordinate array of shape ``[number_of_particles,3]``
        box : float (default 1000)
            Size of simulation in Mpc/*h*.
            
        Returns
        -------
        boolean 
            Array of booleans giving whether each particle is within the shape.
        """
        inside = np.all((pos >= self.posmin) | (pos+box <= self.posmax),axis=1) & np.all((pos <= self.posmax) | (pos-box >= self.posmin),axis=1)
        
        return inside


class Cone():
    """
    Represents a cone. Can check for inclusion of a point and serialize itself in the 
    format accepted by the Spatial3D library's parser.
    
    Parameters
    ----------
    (x, y, z) : float
        The coordinates of the vertex (origin) of the cone in simulation units.
    (dx, dy, dz) : float
        The direction vector; need not be normalized.
    angle : float
        The (half) opening angle of the cone in radians.
    depth : float
        The depth of the cone in simulation units.
    
    """
    def __init__(self,x,y,z,dx,dy,dz,angle,depth):
        
        if (depth <= 0):
            raise ValueError("The depth of a cone must be greater than 0.")
        if ((angle <= 0) or (angle > np.pi)):
            raise ValueError("The angle must be > 0 and <= pi.")
        
        self.vx = x
        self.vy = y
        self.vz = z
        dnorm = np.sqrt(dx*dx+dy*dy+dz*dz)
        self.dx = dx/dnorm
        self.dy = dy/dnorm
        self.dz = dz/dnorm
        self.angle = angle
        self.depth = depth
        self.vert = np.asarray([self.vx,self.vy,self.vz])
        self.dir = np.asarray([self.dx,self.dy,self.dz])
    
        # Could take angle as degrees and convert to radians...
    
    def __str__(self):
        return f"CONE[{self.vx:f},{self.vy:f},{self.vz:f},{self.dx:f},{self.dy:f},{self.dz:f},{self.angle:f},{self.depth:f}]"
    
    def contained(self,pos,box=1000):
        """
        Return a boolean index that will select particles in the shape.

        Parameters
        ----------
        pos : ndarray 
            Particle coordinate array of shape ``[number_of_particles,3]``
        box : float (default 1000)
            Size of simulation in Mpc/*h*.
            
        Returns
        -------
        boolean 
            Array of booleans giving whether each particle is within the shape.
        """
        hbox = box/2
        pdir = pos-self.vert
        pdir[pdir>hbox] -= box
        pdir[pdir<-hbox] += box
        n2 = np.sum(pdir**2,axis=1)
        n = np.sqrt(n2)
        dmag = np.sqrt(np.sum(self.dir**2)) # = 1
        posangle = np.dot(pdir,self.dir)/n/dmag # cosine of angle bet. particles and direction of cone
        inside = (posangle >= np.cos(self.angle)) & (n <= self.depth)
        
        return inside


class ConeSegment():
    """
    Represents a cone segment. Can check for inclusion of a point and serialize itself in the 
    format accepted by the Spatial3D library's parser.
        
    Parameters
    ----------
    (x, y, z) : float
        The coordinates of the vertex (origin) of the cone encompassing this segment in simulation units.
    (dx, dy, dz) : float
        The direction vector, need not be normalized.
    angle : float
        The (half) opening angle of the cone in radians.
    (depthmin, depthmax) : float
        The minimum and maximum depth of the cone segment in simulation units.

    """
    def __init__(self,x,y,z,dx,dy,dz,angle,depthmin,depthmax):
        
        if ((depthmin <= 0) or (depthmin > depthmax)):
            raise ValueError("The min depth of cone segment must be greater than 0 and less than max depth.")
        if ((angle <= 0) or (angle > np.pi)):
            raise ValueError("The angle must be > 0 and <= pi.")

        self.vx = x
        self.vy = y
        self.vz = z
        dnorm = np.sqrt(dx*dx+dy*dy+dz*dz)
        self.dx = dx/dnorm
        self.dy = dy/dnorm
        self.dz = dz/dnorm
        self.angle = angle
        self.depthmin = depthmin
        self.depthmax = depthmax
        self.vert = np.asarray([self.vx,self.vy,self.vz])
        self.dir = np.asarray([self.dx,self.dy,self.dz])
    
    
    def __str__(self):
        return f"CONESEGMENT[{self.vx:f},{self.vy:f},{self.vz:f},{self.dx:f},{self.dy:f},{self.dz:f},{self.angle:f},{self.depthmin:f},{self.depthmax:f}]"
    
    def contained(self,pos,box=1000):
        """
        Return a boolean index that will select particles in the shape.

        Parameters
        ----------
        pos : ndarray 
            Particle coordinate array of shape ``[number_of_particles,3]``
        box : float (default 1000)
            Size of simulation in Mpc/*h*.
            
        Returns
        -------
        boolean 
            Array of booleans giving whether each particle is within the shape.
        """
        hbox = box/2
        pdir = pos-self.vert
        pdir[pdir>hbox] -= box
        pdir[pdir<-hbox] += box
        n2 = np.sum(pdir**2,axis=1)
        n = np.sqrt(n2)
        dmag = np.sqrt(np.sum(self.dir**2)) # = 1
        posangle = np.dot(pdir,self.dir)/n/dmag # cosine of angle bet. particles and direction of cone
        inside = (posangle >= np.cos(self.angle)) & (n <= self.depthmax) & (n >= self.depthmin) 
        
        return inside
