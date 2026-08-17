import numpy as np

def FEAInterpolation (o: np.ndarray, rxyz: np.ndarray, es: np.ndarray, mask: np.ndarray, arr: np.ndarray) -> np.ndarray:
    from Terpdose import _tdfort
    
    loc_o    = np.asfortranarray (o,      dtype=int       )
    loc_rxyz = np.asfortranarray (rxyz,   dtype=np.float64) # Comes into this function as [3, n] because it's made in Fortran
    loc_es   = np.asfortranarray (es + 1, dtype=int       ) # Adjust these for Fortran indexing
    loc_mask = np.asfortranarray (mask,   dtype=bool      )
    loc_arr  = np.asfortranarray (arr,    dtype=np.float64)
    
    return _tdfort.feainterpolation (loc_o, loc_rxyz, loc_es, loc_mask, loc_arr)