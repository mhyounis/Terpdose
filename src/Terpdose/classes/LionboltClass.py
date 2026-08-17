import h5py
import numpy as np
from Terpdose.schema.LionboltSchema import *
from Terpdose.postprocessing.Deposition import *
from Terpdose.classes.PhysicsClasses import XSLibrary
from Terpdose.classes.ParticleClass import Particle
from Terpdose.classes.PhaseSpace import Mesh
from Terpdose.utils.genutils import VersionError, isinstancelist, lbsorted

#  =================
#    LionboltClass  
#  =================

#  ===========================================================================================
#    All indices input and output follow python zero-based indexing. Lionbolt                 
#    HDF5 data is written following Fortran one-based indexing, so the indices                
#    are shifted manually in Terpdose in some cases.                                          
#                                                                                             
#    Index naming conventions: 'energy dof'   --- 'g'  = The energy group                     
#      'angular dof'  --- 'i'  = The discrete ordinate (SN) or real spherical                 
#      harmonic (PN)                                                                          
#      'element'      --- 'e'  = The finite element 'element node' --- 'k'  =                 
#      The local node (a node's index within an element) 'global node'  --- 'kg'              
#      = The global node 'spatial dof'  --- 'sd' = The extended node index,                   
#      i.e., every distinctively indexed                                                      
#                                node in a discontinuous finite element formulation, keeping  
#                                nodes that correspond to the same physical point distinct    
#                                as they exist in different elements.                         
#  ===========================================================================================

# Consider using __enter__ and __exit__ functions, which just allow you to put Lionbolt in a with block, with any tasks upon entering and exit. __enter__ (self) can just be return self, __exit__ (self) can just be self.destroy() (consider the exc_type, exc_value, traceback stuff)

class Lionbolt:
    
    """ The Lionbolt class serves to read in Lionbolt data.
    
    Print an object of this class to view all available attributes.
    
    Parameters
    ----------
    fname : str
        File name corresponding to an HDF5 file output by Lionbolt.
    
    Examples
    --------
    Load a Lionbolt HDF5 file, work with a particle, then destroy:
    
    >>> D = Lionbolt('results.h5')
    >>> el = D.electrons
    >>> D.destroy ()
    
    Attributes
    ----------
    h5 : h5py.File
        h5py.File object corresponding to the user's HDF5 file.
    problem_type : str
        Either 'slab' or 'general' depending on the solve described by this HDF5 file.
    mesh : Mesh
        The mesh used in the solve.
    XSLib : XSLibrary
        Cross section library used in the solve.
    <particle name> : Particle
        Particle solved, named according to how it was written in the HDF5 file.
    """
    
    def __init__ (self, fname):
        
        # Set file path and initialize a blank file, to avoid reading in such a large .h5 file.
        self.h5 = h5py.File(fname, 'r')
        
        # Determine problem type
        self.problem_type = self._get_problem_type ()
        
        # Initialize mesh
        self.mesh = Mesh (self.h5[GROUP_MESH])
        
        # Initialize the XS library
        if GROUP_XSLIBRARY in self.h5.keys():
            self.XSLib = XSLibrary (self.h5[GROUP_XSLIBRARY])
        else:
            self.XSLib = None
        
        # Initialize particles based on what's present
        for group in self.h5.keys():
            if ATTR_OBJTYPE in self.h5[group].attrs:
                objtype = self.h5[group].attrs[ATTR_OBJTYPE].decode('utf-8')
                if objtype == IS_PARTICLE:
                    if self.XSLib is None:
                        XS = None
                    else:
                        XS = self.XSLib[group]
                    new_particle = Particle (
                        self.h5[group],
                        self.mesh,
                        XS # Note, this means a Lionbolt object has duplicate info... is that going to be a problem?
                    )
                    setattr(self, str(group), new_particle)
        
    def _get_problem_type (self):
        
        # Determine the problem type
        
        return self.h5.attrs[ATTR_PROBTYPE].decode('utf-8')
    
    def _close (self):
        
        # Close the HDF5 File
        
        self.h5.close()
    
    def _calculate_deposition (self, deptype):
        
        # Calculates the deposition of a given type.
        
        # First check that deptype is valid. Even though this is a private
        # function it can be useful for me.
        if deptype != 'energy' and deptype != 'dose' and deptype != 'charge':
            raise ValueError (f'LionboltClass.py: _calculate_deposition: {deptype} is not a valid deposition type.')
        
        # Loop through the attributes of self to find all particles, and then calculate their contribution
        for attr in dir(self):
            # Get the attribute
            tmp = getattr(self, attr)
            
            # Check if this is a particle
            if isinstance(tmp, Particle):
                # Calculate the deposition due to this particle
                
                try:
                    fl = tmp.fluence ()
                except:
                    try:
                        fl = tmp.uncollided_fluence ()
                    except:
                        raise ValueError (f'LionboltClass.py: _calculate_deposition: Could not load in <{tmp.name}> fluence nor uncollided fluence.')
                
                if deptype == 'energy':
                    XS = tmp.XS.energy_dep
                elif deptype == 'dose':
                    XS = tmp.XS.dose_dep
                elif deptype == 'charge':
                    XS = tmp.XS.charge_dep
                
                dmap0 = Deposition (fl, self.mesh, XS)
                
                try:
                    # Try to iterate the deposition. If this fails it means deposition has not been initialized
                    dmap += dmap0
                except:
                    # Initialize the deposition since iteration failed
                    dmap = dmap0.copy () # dmap0 will be reassigned, so copy it
        
        return dmap
        
    def energy_deposition (self):
        
        """ Calculates or fetches (if present) the energy deposition map.
        
        Returns
        -------
        dmap : np.float64 [:]
            Energy deposition indexed by spatial degree of freedom.
        """
        
        if DATASET_EDEP + PATTERN_DEP in self.h5.keys():
            # The deposition map was calculated in Lionbolt. Read it.
            dmap = self.h5[DATASET_EDEP + PATTERN_DEP][:]
        else:
            # Must be calculated
            dmap = self._calculate_deposition(deptype='energy')
        
        return dmap
    
    def dose_deposition (self):
        
        """ Calculates or fetches (if present) the dose deposition map.
        
        Returns
        -------
        dmap : np.float64 [:]
            Dose deposition indexed by spatial degree of freedom.
        """
        
        if DATASET_DDEP + PATTERN_DEP in self.h5.keys():
            # The deposition map was calculated in Lionbolt. Read it.
            dmap = np.array(self.h5[DATASET_DDEP + PATTERN_DEP])
        else:
            # Must be calculated
            dmap = self._calculate_deposition(deptype='dose')
        
        return dmap
    
    def charge_deposition (self):
        
        """ Calculates or fetches (if present) the charge deposition map.
        
        Returns
        -------
        dmap : np.float64 [:]
            Charge deposition indexed by spatial degree of freedom.
        """
        
        if DATASET_CDEP + PATTERN_DEP in self.h5.keys():
            # The deposition map was calculated in Lionbolt. Read it.
            dmap = self.h5[DATASET_CDEP + PATTERN_DEP][:]
        else:
            # Must be calculated
            dmap = self._calculate_deposition(deptype='charge')
        
        return dmap
    
    def destroy (self):
        
        """ Deallocates all data 
        """
        
        # Clear metadata
        delattr(self, problem_type)
        
        # Clear the particles
        for attr in dir(self):
            if isinstance(attr, Particle):
                attr.destroy ()
        
        # Clear the mesh
        self.mesh.destroy ()
        
        # Clear the XS Library
        if self.XSLib is not None:
            self.XSLib.destroy ()
        
        # Close the HDF5 file
        self._close ()
    
    #  INFORMATION
    
    def __repr__ (self):
        
        # Lists all available data. Access by printing the object.
        
        from Terpdose.utils.genutils import generic_repr
        
        return generic_repr (self)
