import numpy as np

def LeafElementMask (rg: np.ndarray, c: np.ndarray, o: np.ndarray, es: np.ndarray, xb: np.ndarray, yb: np.ndarray, zb: np.ndarray) -> np.ndarray:
    from Terpdose import _tdfort
    
    loc_rg = np.asfortranarray (rg.T,   dtype=np.float64)
    loc_c  = np.asfortranarray (c + 1,  dtype=int       ) # Adjust these for Fortran indexing
    loc_o  = np.asfortranarray (o,      dtype=int       )
    loc_es = np.asfortranarray (es + 1, dtype=int       ) # Adjust these for Fortran indexing
    loc_xb = np.asfortranarray (xb,     dtype=np.float64)
    loc_yb = np.asfortranarray (yb,     dtype=np.float64)
    loc_zb = np.asfortranarray (zb,     dtype=np.float64)
    
    return _tdfort.leafelementmask (loc_rg, loc_c, loc_o, loc_es, loc_xb, loc_yb, loc_zb)