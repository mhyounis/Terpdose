#!/usr/bin/env python

from __future__ import annotations
from Terpdose import *

def main ():
    
    """ # MHY - WIP, not user-friendly nor documented right now.
    """
    
    import numpy as np
    from argparse import ArgumentParser
    import matplotlib
    import matplotlib.pyplot as plt
    
    #  =====================
    #    Read in arguments  
    #  =====================
    
    parser = ArgumentParser()
    
    parser.add_argument('fname', help='The Lionbolt HDF5 file')
    
    deposition_type = parser.add_mutually_exclusive_group(required=True)
    deposition_type.add_argument('-E', '--energy',  help='Plot (fluence) energy spectrum', action='store_true', default=False)
    deposition_type.add_argument('-A', '--angular', help='Plot angular spectrum',          action='store_true', default=False)
    
    args = parser.parse_args ()
    
    # Make assignments
    h5fname = args.fname
    
    #  ========================
    #    Calculate deposition  
    #  ========================
    
    D = Lionbolt (h5fname)
    
    x = D.electrons.energy.gridpoints ()
    y = D.electrons.fluence_spectrum  ()
    
    matplotlib.use('Agg')
    
    # Enable LaTeX
    plt.rcParams.update({
        'text.usetex' : True,
        'font.family' : 'serif',
        'font.serif'  : ['Computer Modern Roman']
    })
    
    fig, ax = plt.subplots(figsize=(7.1, 4))
    
    ax.stairs (y, x)
    
    pretty_fig (fig, ax)
    
    fig.savefig('figure.png', transparent=False, format='png', bbox_inches='tight', dpi=600)

if __name__ == '__main__':
    main ()