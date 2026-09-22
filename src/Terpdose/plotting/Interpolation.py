import numpy as np
from Terpdose.constants                     import *
from Terpdose.classes.GeometricClasses      import *
from Terpdose.classes.GridClasses           import *
from Terpdose.wrappers.fea.FEAInterpolation import *

#  =================
#    INTERPOLATION  
#  =================

def GeoInterpolation (mesh : Mesh, geo : Line | Plane | Box, arr):
    
    """ This function returns an array which is the initial array interpolated onto 
    the grid within the geometry provided.
    
    It is assumed that arr is indexed by spatial dof.
    
    Parameters
    ----------
    mesh : MeshClass object
        The mesh over which this array is defined.
    geo : GeometricClasses object
        Geometry in which you want to interpolate arr.
    arr : np.float64 [:] or np.float64 [:,:]
        Some quantity defined over spatial d.o.f.
    
    Returns
    -------
    arr_m : np.float64 [:] or np.float64 [:,:] or np.float64 [:,:,:]
        A (masked) array giving the interpolated values on the grid.
        Masking is used when the grid is not rectangular.
    """
    
    grid = Grid (mesh, geo)
    
    if arr.ndim != 1 and arr.ndim != 2:
        raise TypeError ("Interpolation.py: GeoInterpolation: Input array 'arr' must be a 1D or 2D numpy array.")
    
    o = mesh.offset ()
    
    # Begin with flattened arr_m. Note, ref_xyz is already flattened and appropriately shaped
    ogshape = np.shape(grid.xyz)
    
    # Maybe have a logical that checks the shape of element map and mask
    # since they can be carried around either flattened or unflattened
    if arr.ndim == 1:
        # Now interpolate array for each grid point
        arr_m = FEAInterpolation (o, 
                grid.ref_xyz, 
                grid.element_map, 
                grid.element_mask, 
                arr
            )
        
        # Now unflatten
        arr_m = arr_m.reshape(ogshape[:-1])
        
        # Now mask (and return)
        mask = grid.element_mask.reshape(ogshape[:-1])
        
        return np.ma.array(arr_m, mask=mask)
        
    elif arr.ndim == 2:
        
        result = []
        
        for iDir in range(arr.shape[1]):
            arr_slice = arr[:,iDir]
            
            arr_m = FEAInterpolation (o, 
                    grid.ref_xyz, 
                    grid.element_map, 
                    grid.element_mask, 
                    arr_slice
                )
            
            arr_m = arr_m.reshape(ogshape[:-1])
            mask = grid.element_mask.reshape(ogshape[:-1])
            
            masked_slice = np.ma.array(arr_m, mask=mask)
            result.append(masked_slice)
        
        return np.ma.stack (result, axis=0)
        
    else:
        raise ValueError ('Interpolation.py: GeoInterpolation: arr has been given with more than two dimensions')
