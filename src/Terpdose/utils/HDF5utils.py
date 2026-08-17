from __future__ import annotations
import numpy as np

#  ==================
#    HDF5 utilities  
#  ==================

#  ==================================================================================
#    A series of functions which make it easier to pull particularly-formatted data  
#    from HDF5 files.
#    The formatting is, of course, according to Lionbolt and NittanyPhysics output.
#  ==================================================================================

def read_space_angle (h5group : h5py.Group, angles='all'):
    
    '''
        MHY LATER - WIP DESCRIPTION
        
    '''
    
    import h5py
    from Terpdose.utils.genutils import isinstancelist, lbsorted
    
    #  ===================================
    #    Prepare to read in angular data  
    #  ===================================
    
    #  ------------------------------------------
    #    Get HDF5 keys for the different angles  
    #  ------------------------------------------
    
    i_keys = lbsorted(h5group.keys())
    
    NI = len(i_keys)
    
    #  ---------------------------------------------
    #    Determine the set of angular dofs to read  
    #  ---------------------------------------------
    
    if angles == 'all':
        i_list = np.arange(NI)
    else:
        i_list = angles
        
        if (not isinstance (i_list, int)) and (not isinstancelist (i_list, int)):
            raise ValueError ('Only integers and integer lists are allowed as input for angles.')
        
        if isinstance (i_list, int):
            i_list = [i_list]
        
        if (any(i_list) > NI - 1):
            raise ValueError ('An input for angles exceeds the number of angular dofs.')
    
    i_keys = [i_keys[i] for i in i_list]
    
    #  ================================
    #    Begin constructing the array  
    #  ================================
    
    arr = []
    
    for i_key in i_keys:
        arr.append(h5group[i_key])
    
    return np.array(arr)

def ReadPolySpaceAngle (h5 : h5py.Group, energies='all', angles='all'):
    
    '''
        MHY LATER - WIP DESCRIPTION
        
    '''
    
    from Terpdose.utils.genutils import isinstancelist, lbsorted
    
    #  ==================================
    #    Prepare to read in energy data  
    #  ==================================
    
    #  --------------------------------------------
    #    Get HDF5 keys for the different energies  
    #  --------------------------------------------
    
    g_keys = lbsorted(h5.keys())
    
    G = len(g_keys)
    
    #  --------------------------------------------
    #    Determine the set of energy dofs to read  
    #  --------------------------------------------
    
    if energies != 'all':
        g_list = energies
        
        if (not isinstance (g_list, int)) and (not isinstancelist (g_list, int)):
            raise ValueError ('Only integers and integer lists are allowed as input for energies.')
        
        if isinstance (g_list, int):
            g_list = [g_list]
        
        if (any(g_list) > G - 1):
            raise ValueError ('An input for energies exceeds the number of angular dofs.')
        
    else:
        g_list = list(range(0, G))
    
    g_keys = g_keys[g_list]
    
    #  ================================
    #    Begin constructing the array  
    #  ================================
    
    arr = []
    
    for g_str in g_keys:
        
        arr.append(read_space_angle(h5, angles=angles))
    
    return arr