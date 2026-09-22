import numpy as np

def CalculateCurrent (w: np.ndarray, khat: np.ndarray, psi: np.ndarray) -> np.ndarray:
    from Terpdose import _tdfort
    
    loc_w    = np.asfortranarray(w,     dtype=np.float64)
    loc_khat = np.asfortranarray(khat,  dtype=np.float64) # Does not need to be transposed, that is actually more natural within fortran
    loc_psi  = np.asfortranarray(psi.T, dtype=np.float64)
    
    return _tdfort.calculatecurrent (loc_w, loc_khat, loc_psi) # Does not need to be transposed, that is actually more natural within fortran
