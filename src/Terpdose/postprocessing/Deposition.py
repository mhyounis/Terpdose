from __future__ import annotations

def Deposition (fl, mesh : Mesh, depXS):
    
    import numpy as np
    from Terpdose.classes.PhaseSpace import Mesh
    from Terpdose.wrappers.postprocessing.TallyDeposition import TallyDeposition
    
    ''' Calculates a deposition map
        
        (In)  :: fl    [g][sdof] --- Fluence for energy g and spatial dof sdof
        (In)  :: mesh            --- The mesh on which fl exists
        (In)  :: depXS [m, g]    --- The deposition cross section for energy g and material m
        
        (Out) :: arr   [sdof]    --- The deposition at sdof
    '''
    
    arr = np.zeros(mesh.num_spatial_dofs)
    
    # Pre-load mat2sd here outside of the loop
    mat2sd = mesh.mat2sd ()
    
    for (m, g), d in np.ndenumerate(depXS):
        TallyDeposition (mat2sd[m], d, fl[g], arr)
    
    return arr