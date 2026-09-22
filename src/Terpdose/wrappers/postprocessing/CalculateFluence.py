import numpy as np

def CalculateFluence (w: np.ndarray, psi: np.ndarray) -> np.ndarray:
    from Terpdose import _tdfort
    
    loc_w   = np.asfortranarray(w,     dtype=np.float64)
    loc_psi = np.asfortranarray(psi.T, dtype=np.float64) # SLOW? Because you're creating a new array. Maybe just transpose s outside? Not sure what is best here.
    
    return _tdfort.calculatefluence (loc_w, loc_psi)