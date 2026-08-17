import numpy as np

def CalculateCurrent () -> np.ndarray:
    from Terpdose import _tdfort
    
    _tdfort.calculatecurrent ()
