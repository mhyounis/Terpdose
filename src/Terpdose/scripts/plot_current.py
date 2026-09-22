#!/usr/bin/env python

from __future__ import annotations
from Terpdose import *

# WILL LATER SPLIT THIS UP INTO PARSE, MAKE FIG AND AX, AND SAVE FIG. This way making fig and ax can be used as part of Terpdose.

def main ():
    
    """ Plot particle current data from Lionbolt HDF5 files.
    
    To make command-line-interface user-friendly, the plot and geometry options are very few.
    However, you are recommended to use this script as a guideline for making more customizable
    plots.
    
    Can currently only calculate in a 2D plane.
    
    Examples
    --------
    Plot particle named 'electrons', energy group 15 (indexing starts from 1), in YZ plane:
        
        plot_deposition.py results.h5 -p electrons -g 15 --plane YZ
    
    Plot particle named 'photons', sum all energy groups, in XY plane:
    
        plot_deposition.py results.h5 -p photons --plane XY
    
    """
    
    import numpy as np
    from argparse import ArgumentParser
    
    #  =====================
    #    Read in arguments  
    #  =====================
    
    parser = ArgumentParser()
    
    parser.add_argument('fname', help='The Lionbolt HDF5 file')
    
    parser.add_argument('-p', '--particle',
                        help='The name of the particle to be plotted.',
                        required=True,
                        type=str)
    parser.add_argument('-g', '--group', 
                        help="The energy group to be plotted (indexing starts from 1, which has the highest energies). Use 'all' to sum over all groups (default).", 
                        default='all')
    # parser.add_argument('--line-start', 
    #                     help='Plot along a line (for non-slab calculations). Cannot be given with plane nor volume but requires line-end.', 
    #                     nargs=3,
    #                     type=float,
    #                     metavar=('x0', 'y0', 'z0'))
    # parser.add_argument('--line-end', 
    #                     help='Plot along a line (for non-slab calculations). Cannot be given with plane nor volume but requires line-start.', 
    #                     nargs=3,
    #                     type=float,
    #                     metavar=('x1', 'y1', 'z1'))
    parser.add_argument('--plane', 
                        help='Plot along a plane (for non-slab calculations). Cannot be given with line-start, line-end, nor volume.', 
                        choices=['XY', 'XZ', 'YZ'], 
                        default=None)
    # parser.add_argument('--volume', 
    #                     help='Plot a transparency map in full 3D space. Cannot be given with line-start, line-end, nor plane.', 
    #                     action='store_true')
    
    args = parser.parse_args ()
    
    # Make assignments
    h5fname   = args.fname
    p         = args.particle.lower()
    g         = args.group
    # x0        = np.array(args.line_start)
    # x1        = np.array(args.line_end)
    plane     = args.plane
    # do_volume = args.volume
    
    # if do_volume:
    #     try:
    #         import pyvista
    #     except:
    #         raise ModuleNotFoundError ('plot_deposition.py : --volume specified but user does not have pyvista installed.')
    
    # # Determine if there are conflicts
    # if do_volume:
    #     do_line  = False
    #     do_plane = False
    #     
    #     if args.line_start is not None or args.line_end is not None or args.plane is not None:
    #         parser.error ('--volume specified but user also specified either --line-start, --line-end, or --plane.')
    # else:
    #     do_line = args.line_start is not None or args.line_end is not None
    #     
    #     if do_line:
    #         if args.line_start is None or args.line_end is None:
    #             parser.error ('--line-start and --line-end must be specified together.')
    #         if args.plane is not None:
    #             parser.error ('--line-start and --line-end cannot be used with --plane.')
    #     
    #     do_plane = not do_line
    
    #  =====================
    #    Calculate current  
    #  =====================
    
    D = Lionbolt (h5fname)
    
    # Set the problem type
    problem_type = D.problem_type
    
    # if do_volume and problem_type == 'slab':
    #     parser.error ('--volume specified but this .h5 file describes a slab solve.')
    
    # Get the quantity requested
    try:
        particle = getattr(D, p)
    except:
        print(f'Particle {p} not found in Lionbolt data file {h5fname}.')
    
    if g == 'all':
        try:
            # If all groups requested, the current is integrated over all energies (summed over all energy groups)
            J = particle.current ()
            J = sum(J)
        except:
            print(f'Current for particle {p} not able to be calculated from Lionbolt data file {h5fname}.')
    else:
        # try:
        # If only one group requested, the current is taken as the average over the energy group
        J = particle.current (energies=int(g) - 1) # Use fortran indexing for user input
        E = particle.energy.gridpoints ()
        J = J / (E[int(g) - 1] - E[int(g)])
        # except:
        #     print(f'Current for particle {p} not able to be calculated from Lionbolt data file {h5fname}.')
    
    # Plan for the plots
    saveloc =  'current.png'
    title   = r'Particle current'
    if g == 'all':
        # If all groups requested
        cblabel = r'Total particle current per fluence'
    else:
        cblabel = r'Particle current per fluence (MeV$^{-1}$)'
    
    #  =========================
    #    Determine how to plot  
    #  =========================
    
    colors = pretty_patties ()
    
    if problem_type == 'general':
        
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
        fig, ax = plot_2D_vectors (D.mesh, geo, J,
                                   title   = title,
                                   xlabel  = xlabel,
                                   ylabel  = ylabel,
                                   cblabel = cblabel)
        
        # if do_plane:
        #     
        #     nodes  = D.mesh.nodes()
        #     origin = [ HALF * (np.max(nodes[:,0]) - np.min(nodes[:,0])) + np.min(nodes[:,0]),
        #                HALF * (np.max(nodes[:,1]) - np.min(nodes[:,1])) + np.min(nodes[:,1]),
        #                HALF * (np.max(nodes[:,2]) - np.min(nodes[:,2])) + np.min(nodes[:,2]) ]
        #     
        #     if plane == 'XY':
        #         xlabel = r'$x$ (cm)'
        #         ylabel = r'$y$ (cm)'
        #         
        #         sides  = [ np.max(nodes[:,0]) - np.min(nodes[:,0]),
        #                    np.max(nodes[:,1]) - np.min(nodes[:,1]) ]
        #         ax1 = np.array([1.0, 0.0, 0.0]) * sides[0] / TWO
        #         ax2 = np.array([0.0, 1.0, 0.0]) * sides[1] / TWO
        #         
        #     elif plane == 'XZ':
        #         xlabel = r'$x$ (cm)'
        #         ylabel = r'$z$ (cm)'
        #         
        #         sides  = [ np.max(nodes[:,0]) - np.min(nodes[:,0]),
        #                    np.max(nodes[:,2]) - np.min(nodes[:,2]) ]
        #         ax1 = np.array([1.0, 0.0, 0.0]) * sides[0] / TWO
        #         ax2 = np.array([0.0, 0.0, 1.0]) * sides[1] / TWO
        #         
        #     elif plane == 'YZ':
        #         xlabel = r'$y$ (cm)'
        #         ylabel = r'$z$ (cm)'
        #         
        #         sides  = [ np.max(nodes[:,1]) - np.min(nodes[:,1]),
        #                    np.max(nodes[:,2]) - np.min(nodes[:,2]) ]
        #         ax1 = np.array([0.0, 1.0, 0.0]) * sides[0] / TWO
        #         ax2 = np.array([0.0, 0.0, 1.0]) * sides[1] / TWO
        #         
        #     # Now use Terpdose to create the geometry and then the plot
        #     geo = Plane ( n=[500, 500], origin=origin, ax1=ax1, ax2=ax2 )
        #     fig, ax = plot_2D_vectors (D.mesh, geo, J,
        #                                title   = title,
        #                                xlabel  = xlabel,
        #                                ylabel  = ylabel,
        #                                cblabel = cblabel)
            
        # elif do_line:
        #     
        #     geo = Line ( n=100, x0=x0, x1=x1 )
        #     fig, ax = plot_1D (D.mesh, geo, dmap,
        #                        title  = title, 
        #                        xlabel = r'Depth (cm)', 
        #                        ylabel = cblabel,
        #                        color  = colors[0],
        #                        pretty = True)
        # elif do_volume:
        #     
        #     plot_transparency_mapping_pyvista (D.mesh, dmap, cblabel=cblabel, cmap='hot', savelabel=saveloc, prune=prune)
        
    # elif problem_type == 'slab':
    #     
    #     fig, ax = plot_slab (D.mesh, 
    #                          dmap, 
    #                          FMR    = False,
    #                          title  = title,
    #                          xlabel = r'Depth (cm)', 
    #                          ylabel = cblabel,
    #                          color  = colors[0],
    #                          pretty = True)
    
    # if not do_volume:
    #     fig.savefig(saveloc, transparent=False, format='png', bbox_inches='tight', dpi=600)
    
    fig.savefig(saveloc, transparent=False, format='png', bbox_inches='tight', dpi=600)

if __name__ == '__main__':
    main ()