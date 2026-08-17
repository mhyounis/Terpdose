import numpy as np

def FindElement (rg: np.ndarray, c: np.ndarray, o: np.ndarray, es: np.ndarray, r: np.ndarray) -> int:
    from Terpdose import _tdfort
    
    loc_rg = np.asfortranarray (rg.T,   dtype=np.float64)
    loc_c  = np.asfortranarray (c + 1,  dtype=int       ) # Adjust these for Fortran indexing
    loc_o  = np.asfortranarray (o,      dtype=int       )
    loc_es = np.asfortranarray (es + 1, dtype=int       ) # Adjust these for Fortran indexing
    loc_r  = np.asfortranarray (r,      dtype=np.float64)
    
    return _tdfort.findelement (loc_rg, loc_c, loc_o, loc_es, loc_r) - 1