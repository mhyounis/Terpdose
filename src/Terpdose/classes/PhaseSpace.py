import h5py
import numpy as np
from Terpdose.constants import *
from Terpdose.schema.LionboltSchema import *
from Terpdose.utils.genutils import lbsorted
from Terpdose.wrappers.grid.FindElement import *
from Terpdose.wrappers.grid.LeafElementMask import *
from Terpdose.wrappers.fea.MapToReference import *

#  ========================================================================
#      Phase Space Classes                                                 
#                                                                          
#    This module contains several classes related to discretization in a   
#    Lionbolt or NittanyPhysics solve, as well as a class for space-angle  
#    vectors                                                               
#  ========================================================================

#  ========
#    Mesh  
#  ========

class Mesh:
    
    """ This class deals with mesh-based data used by a Lionbolt solve 
    (and thus needed for post-processing and plotting). 
    
    Print an object of this class to view all available attributes.
    
    Parameters
    ----------
    h5group : h5py.Group
        The HDF5 group corresponding to the mesh.
    
    Attributes
    ----------
    h5 : h5py.Group
        The HDF5 group from which this object is built.
    problem_type : str
        Either 'slab' or 'general' depending on the dimensionality of the mesh.
    num_materials : int
        Number of materials present in this mesh.       
    num_elements : int
        Number of finite elements in this mesh.
    num_spatial_dofs : int
        Number of spatial degrees of freedom in this mesh. Note that this is given 
        by visiting every element and adding up the number of nodes it has, even
        if two nodes in different elements describe the same physical point. This is
        consistent with discontinuous Galerkin methods which treat these nodes as
        distinct.
    num_global_nodes : int
        Number of nodes in the global mesh, i.e., number of unique physical points. 
    """
    
    def __init__ (self, h5group : h5py.Group):
        
        self.h5 = h5group
        
        # Initialize some of the attributes
        self.problem_type     = self.h5.attrs[ATTR_PROBTYPE].decode('utf-8')
        self.num_materials    = self.h5.attrs[ATTR_NUMMATS_M][()]
        self.num_elements     = self.h5.attrs[ATTR_NUMELS][()]
        self.num_spatial_dofs = self.h5.attrs[ATTR_NUMSD][()]
        self.num_global_nodes = self.h5.attrs[ATTR_NUMKG][()]
        self.octree_exists    = False # Octree gets generated when requested
    
    def num_element_nodes (self):
        
        """ Gives the number of nodes in an element.
        
        Returns
        -------
        NK : int [:]
            Array giving the number of nodes in the given element.
        
        Examples
        --------
        Let element e describe a first-order tetrahedral element:
        
        >>> print(NK[e])
        4
        """
        
        o = self.offset()
        
        NE = len(o) - 1
        NK = np.empty((NE), dtype=int)
        
        for e in range(0, NE):
            NK[e] = o[e + 1] - o[e]
        
        return NK
    
    def nodes (self):
        
        """ Gives the global mesh.
        
        Returns
        -------
        rg : np.float64 [:,:]
            Array giving the coordinates of the global mesh. Indexed like [node, direction]
        """
        
        arr = np.array(self.h5[DATASET_R])
        
        return arr
    
    def offset (self):
        
        """ Gives the offset array of the mesh. This is a rank-one array that
        allows one to quickly map between an element and local node index to
        the spatial d.o.f. index, namely by taking the element index and
        giving the spatial d.o.f. reached up to that element. See Examples 
        below for a more concrete demonstration. 
        
        Returns
        -------
        o : int [:]
            Offset array indexed by element.
        
        Examples
        --------
        Suppose you want to know the spatial d.o.f. of the kth node found in
        element e, and assign it to variable s. Use:
        
        >>> s = o[e] + k
        """
        
        arr = self.h5[DATASET_OFFSET][:]
        
        # No subtraction is needed to transform to Python indexing because k starting from 0
        # preserves the definition of offset[e], and the indexing internal to offset is adjusted
        # automatically of course.
        
        return arr
    
    def connectivity (self):
        
        """ Gives the connectivity of the mesh. This maps from spatial d.o.f. 
        to global node index.
        
        Returns
        -------
        c : int [:]
            Connectivity array.
        
        Examples
        --------
        Suppose you want to determine the spatial coordinates of element e 
        local node k. Let rg represent the global mesh and o represent the
        offset array. You can use:
        
        >>> s = o[e] + k # Spatial d.o.f. of element e, local node k
        >>> kg = c[s]    # Global node index of this spatial d.o.f.
        >>> print(rg[kg,:]) 
        [1.02, 0.53, 9.03]
        """
        
        return self.h5[DATASET_CONNECTIVITY][:] - 1 # Subtract by 1 to transform Fortran indexing to Python indexing
    
    def volumes (self):
        
        """ Gives the volume of a chosen element. Needed for finite element analysis.
        
        Returns
        -------
        vol : int [:]
            Volume of a given element.
        """
        
        return self.h5[DATASET_VOLS][:]
    
    def extended_nodes (self):
        
        """ Gives the extended nodes, i.e., the mesh nodes indexed by spatial d.o.f.
        
        Returns
        -------
        r : np.float64 [:,:]
            Array giving coordinates of the mesh in terms of the spatial d.o.f. 
            Indexed like [spatial d.o.f., direction]
        """
        
        rg = self.nodes        ()
        c  = self.connectivity ()
        o  = self.offset       ()
        
        r = []
        
        NE = self.num_elements
        for e in range(NE):
            NK = o[e + 1] - o[e]
            for k in range(NK):
                r.append(rg[c[o[e] + k],:])
        
        r = np.array(r)
        
        return r
    
    def mat2sd (self):
        
        """ Maps from a material to the set of spatial d.o.f.s having this material.
        
        Although, note that material assignments are indeed element based, having
        this array be extended to spatial d.o.f.s is a matter of utility for Terpdose's
        post-processing functionalities.
        
        Returns
        -------
        m2sd : list(int [:])
            The set of spatial d.o.f.s with a given material. Indexed like [material][spatial d.o.f.]
        """
        
        m2sd = []
        
        m_keys = lbsorted(self.h5[GROUP_MAT2SD].keys())
        
        for m_str in m_keys:
            m2sd.append(self.h5[f'{GROUP_MAT2SD}/{m_str}'][:] - 1)
        
        return m2sd
    
    def generate_octree (self, max_num_ele=200):
        
        if self.octree_exists:
            return
        
        self.octree_exists = True
        
        self.octree = Octree (mesh=self, max_num_ele=max_num_ele)
    
    def destroy (self):
        
        """ Deallocates all data 
        """
        
        delattr (self, h5)
        delattr (self, problem_type)
        delattr (self, num_materials)
        delattr (self, num_elements)
        delattr (self, num_spatial_dofs)
        delattr (self, num_global_nodes)
        delattr (self, octree)
        self.octree_exists = False
    
    def __repr__ (self):
        
        # Lists all available data. Access by printing the object.
        
        from Terpdose.utils.genutils import generic_repr
        
        return generic_repr (self)

#  ================
#    AngularSpace  
#  ================

class AngularSpace:
    
    """ This class contains angular discretization information.
    
    Print an object of this class to view all available attributes.
    
    Parameters
    ----------
    h5group : h5py.Group
        The HDF5 group corresponding to the particles whose angular 
        discretization is being described.
    
    Attributes
    ----------
    h5 : h5py.Group
        The HDF5 group from which this object is built.
    num_angular_dofs : int
        The number of angular d.o.f.s in this discretization.
        For SN, this is the total number of discrete ordinates. For PN, this
        is the total number of spherical harmonics in the expansion.
    """
    
    def __init__ (self, h5 : h5py.Group):    
        
        self.h5 = h5
        
        # Initialize some of the attributes
        self.num_angular_dofs = self.h5.attrs[ATTR_NUMANGLES][()]
    
    def weights (self):
        
        """ For SN solves, this gives the quadrature weights.
        
        Returns
        -------
        w : np.float64 [:]
            The ith entry in this array is the quadrature weight of abscissa i.
        """
        
        # LATER ON --- will need to verify that this is SN rather than PN
        
        return np.array(self.h5[DATASET_ANGWEIGHTS])
    
    def abscissae (self):
        
        """ For SN solves, this gives the quadrature abscissae.
        
        Returns
        -------
        x : np.float64 [:,:]
            Indexed like [i, direction] for the ith abscissa.
        """
        
        # LATER ON --- will need to verify that this is SN rather than PN
        
        return np.array(self.h5[DATASET_ABSCISSAE]).T
    
    def destroy (self):
        
        """ Deallocates all data 
        """
        
        delattr (self, num_angular_dofs)
    
    def __repr__ (self):
        
        # Lists all available data. Access by printing the object.
        
        from Terpdose.utils.genutils import generic_repr
        
        return generic_repr (self)

class EnergyGrid:
    
    # May do away with this class... seems to have information
    # that is best stored directly with Particle or CrossSections
    
    """ This class contains energy discretization information.
    
    Print an object of this class to view all available attributes.
    
    Parameters
    ----------
    h5group : h5py.Group
        The HDF5 group corresponding to the particles whose energy 
        discretization is being described.
    
    Attributes
    ----------
    h5 : h5py.Group
        The HDF5 group from which this object is built.
    num_energy_dofs : int
        The number of energy d.o.f.s in this discretization.
    """
    
    def __init__ (self, h5 : h5py.File):    
        
        if DATASET_EGRID not in h5:
            return
        
        self.h5 = h5
        
        # Initialize some of the attributes
        self.num_energy_dofs = self.h5.attrs[ATTR_NUMENERGIES][()]
    
    def gridpoints (self):
        
        """ Energy grid points
        
        Returns
        -------
        E : np.float64 [:]
            Gives the set of energy grid points. Not necessarily indexed by energy d.o.f.,
            for instance in the case of multigroup discretization, the size of E will
            actually be num_energy_dofs + 1. The grid points are however still
            organized in descending order such that energy group g is defined by
            energies from E[g + 1] to E[g]
        """
        
        return np.array(self.h5[DATASET_EGRID])
    
    def destroy (self):
        
        """ Deallocates all data 
        """
        
        delattr (self, num_energy_dofs)
    
    def __repr__ (self):
        
        # Lists all available data. Access by printing the object.
        
        from Terpdose.utils.genutils import generic_repr
        
        return generic_repr (self)

class SpaceAngleVector:
    
    """ A class for space-angle vectors, such as angular fluences, sources, etc. 
        Also points to the Mesh and AngularSpace on which a given object is constructed.
        
        Do not try to re-assign or re-populate a variable of this object before 
        destroying.
        
        Print an object of this class to view all available attributes.
        
        Parameters
        ----------
        mesh : Mesh
            The mesh on which this object will be created.
        angular : AngularSpace
            The angular discretization on which this object will be created.
        
        Attributes
        ----------
        <self[i,s]> : np.float64
            Gives the value of the space-angle array at angular d.o.f. i and 
            spatial d.o.f. s
        mesh : Mesh
            The mesh on which this object was created.
        angular : AngularSpace
            The angular discretization on which this object was created.
        num_angular_dofs : int
            Number of angular d.o.f.s.
        num_spatial_dofs : int
            Number of spatial d.o.f.s.
    """
    
    def __init__ (self, mesh : Mesh, angular : AngularSpace):
        
        self.mesh    = mesh
        self.angular = angular
        
        self.num_angular_dofs = angular.num_angular_dofs
        self.num_spatial_dofs = mesh.num_spatial_dofs
        
        self.data = np.empty((self.num_angular_dofs, self.num_spatial_dofs), dtype=np.float64)
    
    def __getitem__ (self, key):
        return self.data[key]
    
    def __setitem__ (self, key, value):
        self.data[key] = value
    
    def __array__(self):
        return self.data.astype(dtype=np.float64, copy=False)
    
    def populate_from_h5 (self, h5group : h5py.Group, angles='all'):
        
        """ Populates this object using space-angle vector data from an HDF5 
        group formatted as per Lionbolt schema conventions.
        
        Parameters
        ----------
        h5group : h5py.Group
            The HDF5 group corresponding to the desired space-angle vector.
        angles : list(int) or int, optional
            The set of angles desired. By default, all are given.
        """
        
        from Terpdose.utils.HDF5utils import read_space_angle
        
        self.data = read_space_angle (h5group, angles=angles) # It is assumed that RHS has the same shape as LHS. Otherwise a user's data is messed up somehow
        
        # If only one angle is requested return a rank one array
        if isinstance(angles, int):
            self.data = self.data[0,:]
        
    def fluence (self):
        
        """ Gives the fluence corresponding to this space-angle vector
        
        Returns
        -------
        fl : np.float64 [:]
            The fluence indexed by spatial d.o.f.
        
        """
        
        from Terpdose.wrappers.postprocessing.CalculateFluence import CalculateFluence
        
        return CalculateFluence (self.angular.weights(), self.data)
        
        # Note to self, the above is inefficient if looped over g...
        # See the below, the original way this was done when fl was
        # assumed to be a list over g, and the SAV was also over g...
        # Must consider whether to incorporate this stuff and where...
        # it depends on use cases. Elsewhere I have a construct like:
        # fl = [angfl[g].fluence() for g in range(G)]
        # Is that super inefficient?
        # Anyway here's the old efficient code:
        #     # According to stack exchange the most efficient way to do this
        #     # will be to pre-allocate by getting the first calculation done manually
        #     f0  = CalculateFluence (x[0].angular.weights(), x[0].data)
        #     arr = np.empty((len(x),) + f0.shape, dtype=np.float64)
        
        #     arr[0,:] = f0
            
        #     for g, xg in enumerate(x[1:], start=1):
        #         arr[g,:] = CalculateFluence (xg.angular.weights(), xg.data)
            
        #     return arr
    
    def destroy (self):
        
        """ Deallocates all data 
        """
        
        self.mesh.destroy()
        self.angular.destroy()
        delattr (self, num_angular_dofs)
        delattr (self, num_spatial_dofs)
        delattr (self, data)
    
    def __repr__ (self):
        
        # Lists all available data. Access by printing the object.
        
        from Terpdose.utils.genutils import generic_repr
        
        return generic_repr (self)

#  ==================
#    OCTREE CLASSES  
#  ==================

# Not phase space but this avoids circular dependency

class Octree: 
    
    ''' This class defines an octree.                                                    
        
        All of the elements in a mesh are grouped into regions ('leaves') of roughly     
        200 elements. Then, grid points are assigned to particular leaves (which is a    
        process that is very efficient), and so in order to find the element in which    
        a grid point lives (for interpolation), one must only search the elements in     
        the leaf in which the grid point lives.                                       '''
    
    def __init__ (
        self, 
        mesh : Mesh,    # Mesh out of which this tree is made
        max_num_ele=200 # Maximum number of elements in a leaf
    ):
        
        self.mesh          = mesh
        self._nodes        = mesh.nodes ()
        self._connectivity = mesh.connectivity ()
        self._offset       = mesh.offset ()
        
        self.xb = np.array([ np.min(self._nodes[:,0]), np.max(self._nodes[:,0]) ]) # It's generally bad practice to label variables like x y z rather than put them in a list.
        self.yb = np.array([ np.min(self._nodes[:,1]), np.max(self._nodes[:,1]) ]) # But I do this because I only have one place where I need to loop over x y z anyway, as well
        self.zb = np.array([ np.min(self._nodes[:,2]), np.max(self._nodes[:,2]) ]) # as for clearer user-input and my own clarity.
        
        # The 'root' is the ancestor leaf
        self.root = Leaf (self, self.xb, self.yb, self.zb)
        # Manually populate as having all elements initially
        self.root.elements = np.arange(self.mesh.num_elements)
        self.num_elements = self.mesh.num_elements
        
        self.root.check_splitting_criterion (max_num_ele=max_num_ele)
    
    def find_element (self, r):
        
        ''' Finds the finite element in which a point r lies. '''
        
        # USE A PATTERN LIKE THIS, JUST FIGURE OUT HOW TO DESCEND THE LEAF IN A SMART WAY
        # ALSO, note that ONE of the parent leaves will hit. So you basically need to check leaf children eight at a time.
        
        leaf = self.root
        
        # Initial check to see if the provided point is even within the leaf structure at all
        if not leaf.point_is_in_leaf (r):
            return -1
        
        # Until you find a leaf whose children don't exist, you
        # loop through the leaf's children and check for the 
        # leaf containing the point r. When you find a leaf
        # containing r, you set it as the parent and continue.
        # When you reach the leaf whose children don't exist, you
        # ought to have an optimal leaf (i.e., one with <= max_num_ele elements)
        # and here you can call the Fortran routine which checks each element for the point.
        # NOTE - I might need to have mesh idx information in terpdose... LATER. As of release
        # only tetrahedra will be validated anyway.
        while leaf.children is not None:
            parent = leaf
            leaf_found = False
            i = -1
            while not leaf_found:
                i = i + 1
                if i == 8:
                    raise ValueError(f'GridClasses.py: Octree: find_element: Leaf containing provided point not found - {[r[0], r[1], r[2]]}')
                
                leaf = parent.children[i]
                
                leaf_found = leaf.point_is_in_leaf (r)
        
        if leaf.num_elements == 0:
            # Shouldn't raise an error - just throw away the grid point. Where/how to handle this? For now record -1
            return -1
        
        return FindElement (self._nodes, self._connectivity, self._offset, leaf.elements, r)

class Leaf:
    
    ''' This class defines a tree leaf with a given mesh and a given region '''
    
    def __init__ (
        self, 
        poctree : Octree,              # Parent octree.
        xb      : tuple[float, float], # x bounds [min, max]
        yb      : tuple[float, float], # y bounds
        zb      : tuple[float, float]  # z bounds
    ):
        
        self.poctree = poctree
        
        self.xb = xb.copy() # Must write these to memory otherwise they can be changed outside of the leaf
        self.yb = yb.copy()
        self.zb = zb.copy()
        
        self.children = None # Initialize as None
        self.elements = None # Don't populate during initialization, because you'll want to be able to control parent elements
    
    def populate (self, parent_elements = None):
        import numpy.ma as ma
        
        ''' Determine the elements in the defined leaf '''
        
        if parent_elements is not None:
            es = parent_elements
        else:
            es = np.arange(self.poctree.mesh.num_elements)
        
        c  = self.poctree._connectivity
        o  = self.poctree._offset
        rg = self.poctree._nodes
        
        lclmask = LeafElementMask (rg, c, o, es, self.xb, self.yb, self.zb)
        
        leafmask = np.zeros (self.poctree.mesh.num_elements, dtype=bool)
        for i, e in enumerate(es):
            leafmask[e] = lclmask[i]
        
        self.elements = np.where(leafmask)[0]
        self.num_elements = len(self.elements)
        
    def split (self):
        
        ''' This function returns eight new leaflets corresponding to    
            splitting the current leaf into eighths.                  '''
        
        # If this leaf has not been initialized, initialize it
        # with all elements present
        if self.elements is None:
            self.populate ()
        
        # Initialize an empty 'leaflets' array
        leaflets = []
        
        # Record half of the total leaf length/width/height
        dx = HALF * (self.xb[1] - self.xb[0])
        dy = HALF * (self.yb[1] - self.yb[0])
        dz = HALF * (self.zb[1] - self.zb[0])
        
        # Record (-, -, -) octant corner of the leaf to serve as an origin
        r0 = np.array([self.xb[0], self.yb[0], self.zb[0]])
        
        for ip in range(8):
            # Un-flatten ip.
            # ip indexes the new leaves, here we map them to grid indices [i,j,k]
            i = ip % 2
            j = (ip // 2) % 2
            k = ip // 4
            
            # Determine the center of the new leaf
            r = np.array([(i + HALF) * dx, (j + HALF) * dy, (k + HALF) * dz]) + r0
            
            # Obtain the new bounds by displacing from the center
            newxb = np.array( [r[0] - HALF * dx, r[0] + HALF * dx] )
            newyb = np.array( [r[1] - HALF * dy, r[1] + HALF * dy] )
            newzb = np.array( [r[2] - HALF * dz, r[2] + HALF * dz] )
            
            # Create the new leaf
            newleaf = Leaf (self.poctree, newxb, newyb, newzb)
            
            # Send in the current leaf's elements as the parent elements for these new leaves
            newleaf.populate (parent_elements=self.elements)
            
            # Append
            leaflets.append(newleaf)
            
        self.children = leaflets
        
    def check_splitting_criterion (self, max_num_ele : int, parent=None):
        
        ''' This function recursively splits leaf children until they have no more than max_num_ele elements. '''
        
        # DO THIS : self.root.check_for_split(max_num_ele)
        #           What this will do is check if self.root.elements is of the right size.
        #           If not, it will itself call self.split(max_num_ele) (self there is root).
        #           Then it will check if the children need to be split via self.check_for_split()
        
        # First see if the elements need to be populated (I'm not sure if there's a case where they don't need to be)
        if self.elements is None:
            # If the parent was NOT provided, this is a root, and was not manually populated elsewhere
            if parent is None:
                self.populate ()
            # Otherwise, use the parent's list of elements to determine the list of elements this leaf has
            # out of a smaller subset
            else:
                self.populate (parent_elements=parent.elements)
        
        self.num_elements = len(self.elements)
        
        if self.num_elements > max_num_ele:
            self.split()
            
            for leaflet in self.children:
                leaflet.check_splitting_criterion (max_num_ele=max_num_ele, parent=self)
        
    def num_lowest_descendents (self):
        
        ''' Returns the number of lowest descendents of the leaf (i.e., descendents with no children). '''
        
        print ('num_lowest_descendents : THIS FUNCTION IS NOT VALIDATED.')
        
        leaf = self
        
        nd = 0
        if leaf.children is not None:
            for i in range(8):
                nd = nd + leaf.children[i].num_lowest_descendents ()
            return nd
        else:
            return 1
    
    def num_total_descendents (self):
        
        ''' Returns the number of total descendents of the leaf. '''
        
        print ('num_total_descendents : THIS FUNCTION IS NOT VALIDATED.')
        
        leaf = self
        
        nd = 1
        if leaf.children is not None:
            for i in range(8):
                nd = nd + leaf.children[i].num_total_descendents ()
            return nd
        else:
            return 1
        
    def point_is_in_leaf (self, r, TOL=1.0e-12):
        
        ''' This function returns True if the point r is in the leaf, False otherwise. '''
        
        return (self.xb[0] - TOL <= r[0] <= self.xb[1] + TOL and
                self.yb[0] - TOL <= r[1] <= self.yb[1] + TOL and
                self.zb[0] - TOL <= r[2] <= self.zb[1] + TOL)