import numpy as np

def MapToReference (rg: np.ndarray, c: np.ndarray, o: np.ndarray, gxyz: np.ndarray, es: np.ndarray, mask: np.ndarray) -> np.ndarray:
    from Terpdose import _tdfort
    
    loc_rg   = np.asfortranarray (rg.T,   dtype=np.float64) # Comes into this function as [NK, 3], needed by Fortran as (3, NK)
    loc_c    = np.asfortranarray (c + 1,  dtype=int       ) # Adjust these for Fortran indexing
    loc_o    = np.asfortranarray (o,      dtype=int       )
    loc_gxyz = np.asfortranarray (gxyz.T, dtype=np.float64) # Comes into this function as [n, 3], needed by Fortran as (3, n)
    loc_es   = np.asfortranarray (es + 1, dtype=int       ) # Adjust these for Fortran indexing
    loc_mask = np.asfortranarray (mask,   dtype=bool      )
    
    return _tdfort.maptoreference (loc_rg, loc_c, loc_o, loc_gxyz, loc_es, loc_mask)