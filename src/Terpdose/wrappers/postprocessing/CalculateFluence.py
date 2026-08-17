import numpy as np

def CalculateFluence (w: np.ndarray, s: np.ndarray) -> np.ndarray:
    from Terpdose import _tdfort
    
    loc_w = np.asfortranarray(w,   dtype=np.float64)
    loc_s = np.asfortranarray(s.T, dtype=np.float64) # SLOW? Because you're creating a new array. Maybe just transpose s outside? Not sure what is best here.
    
    return _tdfort.calculatefluence (loc_w, loc_s)