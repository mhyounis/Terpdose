import numpy as np
from Terpdose.constants import *

# CONSIDER - make these extensions of an abstract class

#  =====================
#    GEOMETRIC CLASSES  
#  =====================

#  ==============================================================
#    These classes define the space in which a grid is created.  
#    That is, while a grid is required to plot, a geometry is    
#    required to create a grid.                                  
#  ==============================================================

# MHY LATER - I need to also create a mesh-slice/line generator. That is to say rather than interpolating onto
# regular grids, I can make it so that I intersect the mesh with, e.g., a plane, and that will generate
# a different 2D mesh on which I can plot much cleaner and more faithfully with something like tricontourf.
# This is definitely something I wanna do soon.

class Line:
    
    """ This class defines and discretizes a 1D region to later be fed to
    create a Grid object.
    
    Parameters
    ----------
    n : int
        Number of bins along the line.
    x0 : np.float64 [3]
        First point along the line.
    x1 : np.float64 [3]
        Final point along the line.
    
    Attributes
    ----------
    xyz : np.float64 [:,:]
        The points on the grid in Cartesian coordinates. Indexed as [gridpoint, direction].
    ehat : np.float64 [3]
        The unit vector pointing from x0 to x1
    """
    
    def __init__ (
        self, 
        n      : int,                        # Number of bins along the line
        x0     : tuple[float, float, float], # First point along the line
        x1     : tuple[float, float, float], # Final point along the line
    ):
        
        if not n > 0:
            raise 'GeometricClasses.py: Line: Provided number of bins is too small. It must be > 0.'
        
        self.n  = n
        self.x0 = np.array (x0)
        self.x1 = np.array (x1)
        
        self.xyz, self.ehat = self._lattice ()
    
    def _lattice (self):
        
        # Internal function to populate grid coordinates.
        
        ehat = (self.x1 - self.x0) / np.linalg.norm(self.x1 - self.x0)
        
        xyz = np.empty((self.n, 3), dtype=np.float64)
        for i in range(0, self.n):
            xyz[i,:] = self.x0 + i * (self.x1 - self.x0) / (self.n - 1)
        
        return xyz, ehat

class Plane:
    
    """ This class defines and discretizes a 2D region to later be fed to
    create a Grid object.
    
    Parameters
    ----------
    n : list(int)
        Number of bins along ax1, ax2.
    origin : np.float64 [3]
        Center of the plane in Cartesian coordinates.
    ax1 : np.float64 [3]
        In-plane axis 1. Magnitude is half-length of plane along axis 1
    ax2 : np.float64 [3]
        In-plane axis 2.
    
    Attributes
    ----------
    xyz : np.float64 [:,:,:]
        The points on the grid in Cartesian coordinates. Indexed as 
        [gridpoint in ax1, gridpoint in ax2, direction].
    bas : np.float64 [:,:]
        The basis vectors of the plane. Indexed as [direction, parametric coordinate]
    """
    
    def __init__ (
        self, 
        n      : tuple[int,   int],          # Number of bins along ax1, ax2.
        origin : tuple[float, float, float], # Center of the plane in Cartesian coordinates.
        ax1    : tuple[float, float, float], # In-plane axis 1. Magnitude is half-length of plane along axis 1
        ax2    : tuple[float, float, float]  # In-plane axis 2.
    ):
        
        if not n[0] > 0:
            raise 'GeometricClasses.py: Plane: Provided number of bins along ax1 is too small. It must be > 0.'
        if not n[1] > 0:
            raise 'GeometricClasses.py: Plane: Provided number of bins along ax2 is too small. It must be > 0.'
        
        self.n      = np.array (n)
        self.origin = np.array (origin)
        self.sides  = np.array ([TWO * np.linalg.norm(np.array(ax1)), TWO * np.linalg.norm(np.array(ax2))])
        self.bas    = np.column_stack((TWO * np.array(ax1) / self.sides[0], TWO * np.array(ax2) / self.sides[1]))
        
        self.xyz = self._lattice ()
    
    def _lattice (self):
        
        # Internal function to populate grid coordinates.
        
        # Alias the parametric coordinate basis
        u = self.bas[:,0]
        v = self.bas[:,1]
        
        du = self.sides[0] / self.n[0]
        dv = self.sides[1] / self.n[1]
        
        # Form the lattice by starting from the bottom left corner (in uv-space) and then incrementing
        xyz0 = self.origin - HALF * (self.sides[0] - du) * u - HALF * (self.sides[1] - dv) * v
        xyz  = np.empty((self.n[0], self.n[1], 3), dtype=np.float64)
        for i in range(0, self.n[0]):
            for j in range(0, self.n[1]):
                xyz[i,j,:] = xyz0 + i * du * u + j * dv * v
        
        return xyz

# class Circle:

# class Shell: # Gives the surface of a sphere. Can be useful for, e.g., GNP studies

class Box:
    
    """ This class defines and discretizes a 3D region to later be fed to
    create a Grid object.
    
    Parameters
    ----------
    n : list(int)
        Number of bins along the ax1, ax2, and ax3.
    origin : np.float64 [3]
        Center of the box in Cartesian coordinates
    ax1 : np.float64 [3]
        In-frame axis 1. Magnitude is half-length of box along axis 1
    ax2 : np.float64 [3]
        In-frame axis 2.
    ax3 : np.float64 [3]
        In-frame axis 3.
    
    Attributes
    ----------
    xyz : np.float64 [:,:,:,:]
        The points on the grid in Cartesian coordinates. Indexed as 
        [gridpoint in ax1, gridpoint in ax2, gridpoint in ax3, direction].
    bas : np.float64 [:,:]
        The basis vectors of the box. Indexed as [direction, parametric coordinate]
    """
    
    def __init__ (
        self, 
        n      : tuple[int,   int,   int],   # Number of bins along the ax1, ax2, and ax3.
        origin : tuple[float, float, float], # Center of the box
        ax1    : tuple[float, float, float], # In-frame axis 1. Magnitude is half-length of box along axis 1
        ax2    : tuple[float, float, float], # In-frame axis 2.
        ax3    : tuple[float, float, float]  # In-frame axis 3.
    ):
        
        if not n[0] > 0:
            raise 'GeometricClasses.py: Box: Provided number of bins along ax1 is too small. It must be > 0.'
        if not n[1] > 0:
            raise 'GeometricClasses.py: Box: Provided number of bins along ax2 is too small. It must be > 0.'
        if not n[2] > 0:
            raise 'GeometricClasses.py: Box: Provided number of bins along ax3 is too small. It must be > 0.'
        
        self.n      = np.array (n)
        self.origin = np.array (origin)
        self.sides  = np.array ([TWO * np.linalg.norm(np.array(ax1)), 
                                 TWO * np.linalg.norm(np.array(ax2)),
                                 TWO * np.linalg.norm(np.array(ax3))])
        self.bas    = np.column_stack((TWO * ax1 / self.sides[0], 
                                       TWO * ax2 / self.sides[1], 
                                       TWO * ax3 / self.sides[2]))
        
        self.xyz = self._lattice ()
    
    def _lattice (self):
        
        # Internal function to populate grid coordinates.
        
        # Set up the parametric coordinate basis
        u = self.bas[:,0]
        v = self.bas[:,1]
        w = self.bas[:,2]
        
        du = self.sides[0] / self.n[0]
        dv = self.sides[1] / self.n[1]
        dw = self.sides[2] / self.n[2]
        
        # Form the lattice by starting from the bottom left corner (in uvw-space) and then incrementing
        xyz0 = self.origin - HALF * (self.sides[0] - du) * u - HALF * (self.sides[1] - dv) * v - HALF * (self.sides[2] - dw) * w
        xyz  = np.empty((self.n[0], self.n[1], self.n[2], 3), dtype=np.float64)
        for i in range(0, self.n[0]):
            for j in range(0, self.n[1]):
                for k in range(0, self.n[2]):
                    xyz[i,j,k,:] = xyz0 + i * du * u + j * dv * v + k * dw * w
        
        return xyz

# class Sphere: