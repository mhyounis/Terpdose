#!/usr/bin/env python

from __future__ import annotations
from Terpdose import *

# WILL LATER SPLIT THIS UP INTO PARSE, MAKE FIG AND AX, AND SAVE FIG. This way making fig and ax can be used as part of Terpdose.

def main ():
    
    """ Plot deposition data from Lionbolt HDF5 files.
    
    To make command-line-interface user-friendly, the plot and geometry options are very few.
    However, you are recommended to use this script as a guideline for making more customizable
    plots.
    
    Command line
    ------------
    User must specify the type of deposition they desire. For non-slab
    calculations, a user must specify either line geometry or plane
    geometry, although plane geometry is limited to XY, YZ, and XZ planes.
    See **Notes** below if you want to use arbitrarily defined planes.
    
    Examples
    --------
    Plot dose in a slab::
    
        plot_deposition.py results.h5 -D
    
    Plot charge deposition along a line from 0 to 10 along the $z-$axis::
    
        plot_deposition.py results.h5 -C --line-start 0.0 0.0 0.0 --line-end 0.0 0.0 10.0
    
    Plot energy deposition in the YZ plane::
    
        plot_deposition.py results.h5 -E --plane YZ
    
    Use ``-h`` or ``--help`` in the command line to list all arguments.
    
    Notes
    -----
    The majority of this script is just handling user input
    to make decisions for the resulting plot. If a user wanted to make
    their own plotting script with more general parameters, this script
    shows how simple it is to use Terpdose to get the deposition map, 
    geometry, and then plot on arbitrary geometries using functions like
    plot_1D / plot_2D / plot_slab.
    """
    
    import numpy as np
    from argparse import ArgumentParser
    
    #  =====================
    #    Read in arguments  
    #  =====================
    
    parser = ArgumentParser()
    
    parser.add_argument('fname', help='The Lionbolt HDF5 file')
    
    deposition_type = parser.add_mutually_exclusive_group(required=True)
    deposition_type.add_argument('-E', '--energy', help='Plot energy deposition', action='store_true', default=False)
    deposition_type.add_argument('-D', '--dose',   help='Plot dose',              action='store_true', default=False)
    deposition_type.add_argument('-C', '--charge', help='Plot charge deposition', action='store_true', default=False)
    
    parser.add_argument('--line-start', 
                        help='Plot along a line (for non-slab calculations). Cannot be given with plane but requires line-end.', 
                        nargs=3,
                        type=float,
                        metavar=('x0', 'y0', 'z0'))
    parser.add_argument('--line-end', 
                        help='Plot along a line (for non-slab calculations). Cannot be given with plane but requires line-start.', 
                        nargs=3,
                        type=float,
                        metavar=('x1', 'y1', 'z1'))
    parser.add_argument('--plane', 
                        help='Plot along a plane (for non-slab calculations). Cannot be given with line-start nor line-end.', 
                        choices=['XY', 'XZ', 'YZ'], 
                        default=None)
    
    args = parser.parse_args ()
    
    # Make assignments
    h5fname     = args.fname
    plot_energy = args.energy
    plot_dose   = args.dose
    plot_charge = args.charge
    x0          = np.array(args.line_start)
    x1          = np.array(args.line_end)
    plane       = args.plane
    
    # Determine if there are conflicts
    do_line = args.line_start is not None or args.line_end is not None
    
    if do_line:
        if args.line_start is None or args.line_end is None:
            parser.error ('--line-start and --line-end must be specified together.')
        if args.plane is not None:
            parser.error ('--line-start and --line-end cannot be used with --plane.')
    
    do_plane = not do_line
    
    #  ========================
    #    Calculate deposition  
    #  ========================
    
    D = Lionbolt (h5fname)
    
    # Set the problem type
    problem_type = D.problem_type
    
    # Get the deposition requested
    if plot_energy:
        saveloc =  'energy_deposition.png'
        title   = r'Energy Deposition'
        cblabel = r'Energy density per fluence (MeV/cm)'
        
        dmap = D.energy_deposition ()
        prune=True
    elif plot_dose:
        saveloc =  'dose_deposition.png'
        title   = r'Dose Deposition'
        cblabel = r'Dose per fluence (MeV-cm$^{2}$/g)'
        
        dmap = D.dose_deposition ()
        prune=True
    elif plot_charge:
        saveloc =  'charge_deposition.png'
        title   = r'Charge Deposition'
        cblabel = r'Charge density per fluence (e/cm)'
        
        dmap = D.charge_deposition ()
        prune=False
    
    #  =========================
    #    Determine how to plot  
    #  =========================
    
    if problem_type == 'general':
        
        if do_plane:
            
            nodes  = D.mesh.nodes()
            origin = [ HALF * (np.max(nodes[:,0]) - np.min(nodes[:,0])) + np.min(nodes[:,0]),
                       HALF * (np.max(nodes[:,1]) - np.min(nodes[:,1])) + np.min(nodes[:,1]),
                       HALF * (np.max(nodes[:,2]) - np.min(nodes[:,2])) + np.min(nodes[:,2]) ]
            
            if plane == 'XY':
                xlabel = r'$x$ (cm)'
                ylabel = r'$y$ (cm)'
                
                sides  = [ np.max(nodes[:,0]) - np.min(nodes[:,0]),
                           np.max(nodes[:,1]) - np.min(nodes[:,1]) ]
                ax1 = np.array([1.0, 0.0, 0.0]) * sides[0] / TWO
                ax2 = np.array([0.0, 1.0, 0.0]) * sides[1] / TWO
                
            elif plane == 'XZ':
                xlabel = r'$x$ (cm)'
                ylabel = r'$z$ (cm)'
                
                sides  = [ np.max(nodes[:,0]) - np.min(nodes[:,0]),
                           np.max(nodes[:,2]) - np.min(nodes[:,2]) ]
                ax1 = np.array([1.0, 0.0, 0.0]) * sides[0] / TWO
                ax2 = np.array([0.0, 0.0, 1.0]) * sides[1] / TWO
                
            elif plane == 'YZ':
                xlabel = r'$y$ (cm)'
                ylabel = r'$z$ (cm)'
                
                sides  = [ np.max(nodes[:,1]) - np.min(nodes[:,1]),
                           np.max(nodes[:,2]) - np.min(nodes[:,2]) ]
                ax1 = np.array([0.0, 1.0, 0.0]) * sides[0] / TWO
                ax2 = np.array([0.0, 0.0, 1.0]) * sides[1] / TWO
                
            # Now use Terpdose to create the geometry and then the plot
            geo = Plane ( n=[500, 500], origin=origin, ax1=ax1, ax2=ax2 )
            fig, ax = plot_2D (D.mesh, geo, dmap,
                               title=title,
                               xlabel=xlabel,
                               ylabel=ylabel,
                               cblabel=cblabel,
                               prune=prune)
        else:
            
            geo = Line ( n=100, x0=x0, x1=x1 )
            fig, ax = plot_1D (D.mesh, geo, dmap,
                               title=title, 
                               xlabel=r'Depth (cm)', 
                               ylabel=cblabel,
                               pretty=True)
        
    elif problem_type == 'slab':
        
        fig, ax = plot_slab (D.mesh, 
                             dmap, 
                             FMR=False,
                             title=title,
                             xlabel=r'Depth (cm)', 
                             ylabel=cblabel,
                             pretty=True)
    
    fig.savefig(saveloc, transparent=False, format='png', bbox_inches='tight', dpi=600)

if __name__ == '__main__':
    main ()