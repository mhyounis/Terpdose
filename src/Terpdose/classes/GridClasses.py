import numpy as np
from Terpdose.constants import *
from Terpdose.classes.PhaseSpace import Mesh
from Terpdose.classes.GeometricClasses import *
from Terpdose.wrappers.grid.FindElement import *
from Terpdose.wrappers.grid.LeafElementMask import *
from Terpdose.wrappers.fea.MapToReference import *

#  ================
#    GRID CLASSES  
#  ================

# MHY LATER - rework. The octree should belong to the mesh. This makes grid obsolete. GOOD.
#             scatter grid's uses in geometry. Then maybe rename GeometricClasses to grid or something.
#             This makes geo object have the element interpolation info. Ultimate point is that grid is a ridiculous
#             object to require in addition to geo.

class Grid:
    
    """ This class is largely internal and may be removed shortly. It is just a connection point between
    the user-friendly geometry classes and mesh class. It used to be more important, storing the octree
    and all interpolation information, but it was decided that the mesh should store the octree, making
    this class contain a bunch of important and necessary procedures that I would like to transfer
    to the geometry object but I'm not sure how I will do that.
    
    A cleanup is imminent but for now try not to use this class directly. It is not needed by the user.
    """
    
    def __init__ (self, mesh : Mesh, geo : Line | Plane | Box, max_num_ele=200):
        
        ''' Create the grid '''
        
        self.mesh = mesh
        
        mesh.generate_octree ()
        
        # offset, connectivity, and global nodes are cached because they are used so often
        self._nodes        = mesh.nodes        ()
        self._connectivity = mesh.connectivity ()
        self._offset       = mesh.offset       ()
        
        self.set_geo (geo)
        
    def _element_map (self, shape='flat'):
        
        ''' Determines the map from grid point to element A cull mask is given
            by element_map == -1. This means, if mask = True, the grid point
            does NOT live within an element (this mask convention follows
            numpy.ma).
            
            Can choose to return the grid-shaped arrays if shape='grid' is
            provided.
        '''
        
        if shape != 'flat' and shape != 'grid':
            raise ValueError('GridClasses.py:_element_map: Unknown shape provided - ' + shape)
        
        # First, flatten xyz, because it may be provided in 1D, 2D, or 3D.
        ogshape  = np.shape(self.xyz)
        self.xyz = self.xyz.reshape(-1, self.xyz.shape[-1])
        
        # Now find the element in which each grid point lies
        emap = np.zeros (self.xyz.shape[0], dtype=int)
        for i, r in enumerate(self.xyz):
            emap[i] = self.mesh.octree.find_element (r)
        # Finally, reshape xyz back to its original shape, and shape emap to it
        # if shape='grid'.
        self.xyz = self.xyz.reshape (ogshape)
        if shape=='grid':
            emap = emap.reshape(ogshape[:-1]) # Removes the last rank from ogshape
        
        emask = emap == -1
        
        return emap, emask
    
    def _get_ref_xyz (self):
        
        ''' Returns the grid, where each point has been mapped to its elements reference element. '''
        
        ''' NOTE LATER --- shape is [3, n] '''
        
        # First flatten grid
        ogshape = np.shape(self.xyz)
        self.xyz = self.xyz.reshape(-1, self.xyz.shape[-1])
        
        # Get the needed arrays.
        rg = self._nodes
        c  = self._connectivity
        o  = self._offset
        
        # Now send to Fortran
        ref_xyz = MapToReference (rg, c, o, self.xyz, self.element_map, self.element_mask)
        
        # Unflatten grid
        self.xyz = self.xyz.reshape (ogshape)
        
        return ref_xyz
    
    def set_geo (self, geo : Line | Plane | Box):
        
        """ Set / change the geometry to be used for this grid. Useful when you don't want to recreate    
            an octree (which takes some time).
        
        Parameters
        ----------
        geo : GeometricClasses object
        """
        
        self.geo = geo
        self.xyz = self.geo.xyz
        
        self.element_map, self.element_mask = self._element_map (shape='flat')
        
        self.ref_xyz = self._get_ref_xyz ()
