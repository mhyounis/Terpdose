import numpy as np

def TallyDeposition (mat2sd: np.ndarray, depXS: float, fl: np.ndarray, DEP: np.ndarray) -> np.ndarray:
    from Terpdose import _tdfort
    
    loc_mat2sd = np.asfortranarray(mat2sd, dtype=np.int32) + 1 # Add 1 for Fortran indexing
    loc_fl     = np.asfortranarray(fl,   dtype=np.float64)
    DEP        = np.asfortranarray(DEP, dtype=np.float64)
    
    _tdfort.tallydeposition (loc_mat2sd, depXS, loc_fl, DEP) # Double check this with the return statement... figure out best pattern for when Fortran is iterating on an array
