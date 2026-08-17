from __future__ import annotations
import numpy as np
from Terpdose.constants import *
from Terpdose.schema.LionboltSchema import *
from Terpdose.utils.genutils import lbsorted

#  ====================
#      Particle Class  
#  ====================

class Particle:
    
    """ The Particle class stores particle-specific data.
    
    Print an object of this class to view all available attributes.
    
    Parameters
    ----------
    h5group : h5py.Group
        The HDF5 group corresponding to the particle.
    mesh : Mesh
        The mesh used in the solve.
    XS : CrossSections
        The cross sections of the given particle
    
    Attributes
    ----------
    h5 : h5py.Group
        The HDF5 group from which this object is built.
    name : str
        The name of the particle according to the HDF5 file.
    XS : CrossSections
        The cross sections of the given particle.
    angular : AngularSpace
        The angular discretization of the given particle.
    energy : EnergyGrid
        The energy discretization of the given particle.
    <DATASET> : arbitrary type, rank
        Any non-standard datasets found within this particle's HDF5 group will be accessible. 
        HDF5 attributes will be ignored however and must be accessed manually by the user.
    """
    
    def __init__ (self, h5group : h5py.Group, mesh : Mesh, XS : CrossSections):
        
        from Terpdose.classes.PhaseSpace import Mesh, AngularSpace, EnergyGrid
        from Terpdose.classes.PhysicsClasses import CrossSections
        
        # Assign the h5 group
        self.h5 = h5group
        
        # MHY LATER - NOW VALIDATE THAT THIS IS A VALID PARTICLE h5 BLOCK.
        # MAYBE USE ATTRIBUTES?
        
        # Assign the particle name
        self.name = h5group.name.split('/')[-1]
        
        # Assign the mesh # I don't like how this remains assigned and is technically visible to users. I also worry about copying
        self._mesh = mesh
        
        # Assign the XS # Check for None?
        self.XS = XS
        
        # Initialize phase space attributes
        self.angular = AngularSpace (self.h5)
        self.energy  = EnergyGrid   (self.h5)
    
    def __getattr__(self, DATASET):
        
        # This will allow a user to fetch an arbitrary dataset
        # by just doing like self.DATASET
        
        try:
            return self._get_dataset(DATASET)
        except KeyError:
            raise AttributeError(DATASET)
    
    def angular_fluence (self, energies='all', angles='all'):
        
        """ Fetches (if present) the angular fluence.
        
        Parameters
        ----------
        energies : list(int) or int, optional
            The set of energy indices for which you want to get the angular fluence.
            By default, all are provided.
        angles : list(int) or int, optional
            The set of angular indices for which you want to get the angular fluence.
            By default, all are provided.
        
        Returns
        -------
        angfl : list(SpaceAngleVector) or SpaceAngleVector or None
            The angular fluence. If multiple energies are requested, then this is a
            list of SpaceAngleVector objects. Otherwise it is a single SpaceAngleVector.
            If the HDF5 file does not contain angular fluences, None is returned.
        """
        
        return self._get_polyenergetic_space_angle_vector (GROUP_ANG_FL,
                                                           energies=energies,
                                                           angles=angles)
    
    def uncollided_angular_fluence (self, energies='all', angles='all'):
        
        """ Fetches (if present) the uncollided angular fluence.
        
        Parameters
        ----------
        energies : list(int) or int, optional
            The set of energy indices for which you want to get the angular fluence.
            By default, all are provided.
        angles : list(int) or int, optional
            The set of angular indices for which you want to get the angular fluence.
            By default, all are provided.
        
        Returns
        -------
        angfl : list(SpaceAngleVector) or SpaceAngleVector or None
            The uncollided angular fluence. If multiple energies are requested, then this is a
            list of SpaceAngleVector objects. Otherwise it is a single SpaceAngleVector.
            If the HDF5 file does not contain uncollided angular fluences, None is returned.
        """
        
        return self._get_polyenergetic_space_angle_vector (GROUP_ANG_FL_UNC,
                                                           energies=energies,
                                                           angles=angles)
    
    def fluence (self, energies='all'):
        
        """ Calculates or fetches (if present) the fluence.
        
        Parameters
        ----------
        energies : list(int) or int, optional
            The set of energy indices for which you want to get the angular fluence.
            By default, all are provided.
        
        Returns
        -------
        fl : list(np.float64 [:]) or np.float64 or None
            The fluence. If multiple energies are requested, then this is a
            list of rank-one np.float64 arrays (indexed by spatial d.o.f.). Otherwise 
            it is a single np.float64 array. If the HDF5 file does not contain fluences
            nor the angular fluences, None is returned.
        """
        
        # First check if fluences are present in the H5 file
        
        if GROUP_FL in self.h5:
            # If it's present, read it
            fl = self._get_polyenergetic_space_vector (GROUP_FL, energies=energies)
        else:
            # If not, calculate it
            # We go group-by-group to avoid reading in the entire angular
            # fluence set.
            fl = []
            for g in energies:
                angfl = self.angular_fluence(energies=g)
                
                fl.append(angfl.fluence())
            
            # Convert to single np array if only one energy was requested
            if len(energies) == 1:
                fl = fl[0]
        
        return fl
    
    def uncollided_fluence (self, energies='all'):
        
        """ Calculates or fetches (if present) the uncollided fluence.
        
        Parameters
        ----------
        energies : list(int) or int, optional
            The set of energy indices for which you want to get the angular fluence.
            By default, all are provided.
        
        Returns
        -------
        fl : list(np.float64 [:]) or np.float64 or None
            The uncollided fluence. If multiple energies are requested, then this is a
            list of rank-one np.float64 arrays (indexed by spatial d.o.f.). Otherwise 
            it is a single np.float64 array. If the HDF5 file does not contain fluences
            nor the angular fluences, None is returned.
        """
        
        # First check if fluences are present in the H5 file
        
        if GROUP_FL_UNC in self.h5:
            # If it's present, read it
            fl = self._get_polyenergetic_space_vector (GROUP_FL_UNC, energies=energies)
        else:
            # If not, calculate it
            # We go group-by-group to avoid reading in the entire angular
            # fluence set.
            fl = []
            for g in energies:
                angfl = self.uncollided_angular_fluence(energies=g)
                
                fl.append(angfl.fluence())
            
            # Convert to single np array if only one energy was requested
            if len(energies) == 1:
                fl = fl[0]
        
        return fl
    
    def _get_polyenergetic_space_angle_vector (self, GROUP, energies='all', angles='all'):
        
        # For more general use, creates a list of space-angle vectors (listed over energy) from the 
        # group named 'GROUP.'
        
        # MHY LATER - Check whether or not GROUP contains the needed data. Raise error or give None?
        # If giving None go back to LionboltClass/_calculate_deposition and change that to check
        # if fl is None.
        
        # Get energy keys
        g_keys = lbsorted(self.h5[GROUP].keys())
        
        # Determine how many energies the user is requesting, and which
        if energies == 'all':
            G      = len(g_keys)
            g_list = np.arange(G)
        else:
            if isinstance(energies, int):
                G = 1
                g_list = [energies]
            elif isinstance(energies, list):
                if not isinstancelist(energies, int):
                    raise TypeError ("'energies' argument must be either 'all', an integer, or a list of integers.")
                G      = len(energies)
                g_list = energies
        
        g_keys = [g_keys[g] for g in g_list]
        
        # Now construct the list of space-angle vectors
        savg = []
        for g in g_keys:
            savg.append(self._get_space_angle_vector(GROUP=f'{GROUP}/{g}', angles=angles))
        
        # Pull out of list if energies is just one integer
        if isinstance(energies, int):
            savg = savg[0]
        
        return savg
    
    def _get_space_angle_vector (self, GROUP, angles='all'):
        
        # For more general use, pulls a space-angle vector from the group named 'GROUP.'
        
        from Terpdose.classes.PhaseSpace import SpaceAngleVector
        
        # MHY LATER - Check whether or not GROUP contains the needed data
        
        sav = SpaceAngleVector(self._mesh, self.angular)
        sav.populate_from_h5(self.h5[GROUP], angles=angles)
        
        return sav
    
    def _get_polyenergetic_space_vector (self, GROUP, energies='all'):
        
        # For more general use, pulls a space vector from the group named 'GROUP.'
        
        # Get energy keys
        g_keys = lbsorted(self.h5[GROUP].keys())
        
        # Determine how many energies the user is requesting, and which
        if energies == 'all':
            G      = len(g_keys)
            g_list = np.arange(G)
        else:
            if isinstance(energies, int):
                G = 1
                g_list = [energies]
            elif isinstance(energies, list):
                if not isinstancelist(energies, int):
                    raise TypeError ("ParticleClass.py: _get_polyenergetic_space_vector: 'energies' argument must be either 'all', an integer, or a list of integers.")
                G      = len(energies)
                g_list = energies
        
        g_keys = [g_keys[g] for g in g_list]
        
        # Now construct the list of space-angle vectors
        savg = []
        for g in g_keys:
            savg.append(np.array(self.h5[f'{GROUP}/{g}']))
        
        # Pull out of list if energies is just one integer
        if isinstance(energies, int):
            savg = savg[0]
        
        return savg
    
    def _get_dataset (self, GROUP):
        
        # Gets a dataset from the HDF5 file.
        
        return self.h5[GROUP][:]
    
    def __repr__ (self):
        
        # Lists all available data. Access by printing the object.
        
        from Terpdose.utils.genutils import generic_repr
        
        return generic_repr (self)