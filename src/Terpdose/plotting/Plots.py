from __future__ import annotations

# POSSIBLE PLAN --- Make these functions accept mesh, array, geo,
# then they figure out how to actually plot.
# Or just make a 1D plot 2D plot etc.
# Could still have wrapper functions that

def pretty_patties ():
    
    """ Gives a color scheme for curves, as a list of hex values.
    """
    
    return [
        '#009cde', # Lionbolt blue
        '#585858', # NittanyPhysics gray
        '#07c100', # OpenRPS green
        '#363636', # NittanyPhysics deep gray
        '#004e6f', # Lionbolt deep blue
        '#004436', # OpenRPS deep green
        '#52ccff', # Lionbolt light blue
    ]

def pretty_fig (fig, ax):
    
    """ Formats 1D figures in a specific way. Users can control whether or not this
    is done by using the optional boolean 'pretty' in the 1D plotting functions in
    Terpdose.
    
    A user may wish to turn this off during Terpdose's plotting and then apply this
    manually after plotting multiple curves.
    
    A user may also wish to tweak some of the linewidths and other sizes. In that case
    they should simply copy this function from the source code and make their own function.
    
    In the future there may be an option which cycles between a set of size parameters,
    with user input being like 'small', 'medium', 'large'.
    """
    
    ax.grid(which='major', color='lightgray', linestyle='-', linewidth=0.15)
    ax.grid(which='minor', color='lightgray', linestyle=':', linewidth=0.1, dashes=(6, 6))
    
    ax.minorticks_on()
    
    for spine in ax.spines.values():
        spine.set_linewidth(0.2)
        spine.set_color('black')
        spine.set_zorder(10)
    
    ax.tick_params(axis='both', which='major', direction='in', color='black', 
                   length=4, width=0.2, top=True, right=True, zorder=0)
    ax.tick_params(axis='both', which='minor', direction='in', color='black', 
                   length=2, width=0.2, top=True, right=True, zorder=0)
    
    return fig, ax

def plot_1D (mesh, geo : Line, array, title=None, xlabel=r'$d$ (cm)', ylabel=None, color='blue', pretty=True, dots=False):
    
    """ Plots on a 1D geometry given by geo.
    
    Parameters
    ----------
    mesh : Mesh
        The mesh in which to solve.
    geo : Line
        The 1D geometry over which to plot.
    array : np.float64 [:]
        Some quantity defined over spatial d.o.f.
    title : str, optional
        Title of the plot. By default, no title is used.
    xlabel : str, optional
        x-label of the plot. By default, r'$d$ (cm)'
    ylabel : str, optional
        y-label of the plot. By default, no label is used.
    pretty : boolean, optional
        Whether or not to format the figures in a particular 'pretty' way. 
        By default, this is done, but if a user wishes to use this function to, e.g.,
        plot multiple curves, they may wish to turn this OFF and then apply it themselves
        at the very end by sending fig and ax to the function pretty_fig (fig, ax).
    dots : boolean, optional
        Whether or not to use dots on the curve.
    
    Returns
    -------
    fig : matplotlib.Figure
        matplotlib figure.
    ax : matplotlib.Axes
        matplotlib axes.
    """
    
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from Terpdose import Mesh, GeoInterpolation
    
    matplotlib.use('Agg')
    
    # Enable LaTeX
    plt.rcParams.update ({
        'text.usetex' : True,
        'font.family' : 'serif',
        'font.serif'  : ['Computer Modern Roman']
    })
    
    u = np.array([np.dot(geo.xyz[i,:] - geo.x0, geo.ehat) for i in range(geo.n)])
    y = GeoInterpolation (mesh, geo, array)
    
    fig, ax = plt.subplots(figsize=(7.1, 4))
    
    if dots:
        line = '-o'
    else:
        line = '-'
    
    ax.plot (u, y, line, markersize=0.5, linewidth=1.0, color=color)
    
    if title is not None:
        ax.set_title(title, fontsize=12, pad=10)
    ax.set_xlabel(xlabel, fontsize=10)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=10)
    ax.set_xlim(u.min(), u.max())
    
    if pretty:
        pretty_fig (fig, ax)
    
    return fig, ax

def plot_2D (mesh, geo : Plane, array, title=None, xlabel=r'$x$ (cm)', ylabel=r'$y$ (cm)', cblabel=None, prune=False):
    
    """ Plots on a 2D geometry given by geo.
    
    Parameters
    ----------
    mesh : Mesh
        The mesh in which to solve.
    geo : Plane
        The 2D geometry over which to plot.
    array : np.float64 [:]
        Some quantity defined over spatial d.o.f.
    title : str, optional
        Title of the plot. By default, no title is used.
    xlabel : str, optional
        x-label of the plot. By default, r'$x$ (cm)'
    ylabel : str, optional
        y-label of the plot. By default, r'$y$ (cm)'
    cblabel : str, optional
        colorbar label of the plot. By default, no label is used.
    prune : boolean, optional
        Whether or not to take negative values to zero. Further notes on this
        option are in the source code. By default, this is false.
    
    Returns
    -------
    fig : matplotlib.Figure
        matplotlib figure.
    ax : matplotlib.Axes
        matplotlib axes.
    """
    
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from Terpdose import Mesh, GeoInterpolation
    
    matplotlib.use('Agg')
    
    # Enable LaTeX
    plt.rcParams.update ({
        'text.usetex' : True,
        'font.family' : 'serif',
        'font.serif'  : ['Computer Modern Roman']
    })
    
    U = np.sum((geo.xyz - geo.origin) * geo.bas[:,0], axis=-1)
    V = np.sum((geo.xyz - geo.origin) * geo.bas[:,1], axis=-1)
    
    extent = (U.min(), U.max(), V.min(), V.max())
    
    z = GeoInterpolation (mesh, geo, array)
    
    # Get rid of negative values if requested.
    # In principle, and as we have found during validation, the negative values
    # should be extremely few and very small in magnitude for physical
    # quantities that do not expect them (energy, dose, etc.).
    # Still, even just one negative value can throw off the color bar
    # and thus make the plot harder to read.
    # Take care not to turn this on when plotting something like charge.
    if prune:
        z[z < 0] = 0.0
    
    fig, ax = plt.subplots(figsize=(7.1, 4))
    
    im = ax.imshow (
        z.T,
        origin='lower',
        extent=extent,
        cmap='jet',
        interpolation='none',
        aspect='equal',
    )
    
    if cblabel is None:
        fig.colorbar (im, ax=ax)
    else:
        fig.colorbar (im, ax=ax, label=cblabel)
    
    if title is not None:
        ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=14)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=14)
    
    return fig, ax

def plot_2D_vectors (mesh, geo : Plane, array, title=None, xlabel=r'$x$ (cm)', ylabel=r'$y$ (cm)', cblabel=None):
    
    """ Plots vectors on a 2D geometry given by geo.
    
    Parameters
    ----------
    mesh : Mesh
        The mesh in which to solve.
    geo : Plane
        The 2D geometry over which to plot.
    array : np.float64 [:,:]
        Set of vectors, indexed like [spatial d.o.f., Cartesian component]
    title : str, optional
        Title of the plot. By default, no title is used.
    xlabel : str, optional
        x-label of the plot. By default, r'$x$ (cm)'
    ylabel : str, optional
        y-label of the plot. By default, r'$y$ (cm)'
    cblabel : str, optional
        colorbar label of the plot. By default, no label is used.
    
    Returns
    -------
    fig : matplotlib.Figure
        matplotlib figure.
    ax : matplotlib.Axes
        matplotlib axes.
    """
    
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from Terpdose import Mesh, GeoInterpolation
    
    matplotlib.use('Agg')
    
    # Enable LaTeX
    plt.rcParams.update ({
        'text.usetex' : True,
        'font.family' : 'serif',
        'font.serif'  : ['Computer Modern Roman']
    })
    
    U = np.sum((geo.xyz - geo.origin) * geo.bas[:,0], axis=-1)
    V = np.sum((geo.xyz - geo.origin) * geo.bas[:,1], axis=-1)
    
    extent = (U.min(), U.max(), V.min(), V.max())
    
    vecs = GeoInterpolation (mesh, geo, array)
    
    vmap = np.linalg.norm(vecs, axis=0)
    
    v1 = np.einsum('ijk,i->jk', vecs, geo.bas[:,0]) / vmap
    v2 = np.einsum('ijk,i->jk', vecs, geo.bas[:,1]) / vmap
    
    fig, ax = plt.subplots(figsize=(7.1, 4))
    
    im = ax.imshow (
        vmap.T,
        origin='lower',
        extent=extent,
        cmap='jet',
        interpolation='none',
        aspect='equal',
    )
    
    # DECIMATING THE GRID FOR QUIVER... you don't want a mess of arrows
    skip = 25 # Number of vectors to skip... should get a good sense of what is generally best. Probably use size of slice to determine and fix the number of vectors along each dimension instead
    
    x = np.linspace(extent[0], extent[1], vmap.shape[0])
    y = np.linspace(extent[2], extent[3], vmap.shape[1])
    
    X, Y = np.meshgrid(x, y)
    
    ax.quiver (
        X   [::skip, ::skip],
        Y   [::skip, ::skip],
        v1.T[::skip, ::skip],
        v2.T[::skip, ::skip],
        color='white',
        pivot='mid'
    )
    
    if cblabel is None:
        fig.colorbar (im, ax=ax)
    else:
        fig.colorbar (im, ax=ax, label=cblabel)
    
    if title is not None:
        ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel, fontsize=14)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=14)
    
    return fig, ax

def plot_slab (mesh : Mesh, array, FMR=False, title=None, xlabel=r'Depth (cm)', ylabel=None, color='blue', pretty=True, dots=False):
    
    """ Plots on a full slab mesh.
    
    Parameters
    ----------
    mesh : Mesh
        The mesh in which to solve.
    geo : Line
        The 1D geometry over which to plot.
    array : np.float64 [:]
        Some quantity defined over spatial d.o.f.
    title : str, optional
        Title of the plot. By default, no title is used.
    xlabel : str, optional
        x-label of the plot. By default, r'$d$ (cm)'
    ylabel : str, optional
        y-label of the plot. By default, no label is used.
    pretty : boolean, optional
        Whether or not to format the figures in a particular 'pretty' way. 
        By default, this is done, but if a user wishes to use this function to, e.g.,
        plot multiple curves, they may wish to turn this OFF and then apply it themselves
        at the very end by sending fig and ax to the function pretty_fig (fig, ax).
    dots : boolean, optional
        Whether or not to use dots on the curve.
    
    Returns
    -------
    fig : matplotlib.Figure
        matplotlib figure.
    ax : matplotlib.Axes
        matplotlib axes.
    """
    
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from Terpdose import Mesh
    
    matplotlib.use('Agg')
    
    # Enable LaTeX
    plt.rcParams.update ({
        'text.usetex' : True,
        'font.family' : 'serif',
        'font.serif'  : ['Computer Modern Roman']
    })
    
    z = mesh.extended_nodes()[:,0]
    
    if FMR:
        z = z / max(z)
    
    fig, ax = plt.subplots(figsize=(7.1, 4))
    
    if dots:
        line = '-o'
    else:
        line = '-'
    
    ax.plot (z, array, line, markersize=0.5, linewidth=1.0, color=color)
    
    if title is not None:
        ax.set_title(title, fontsize=12, pad=10)
    ax.set_xlabel(xlabel, fontsize=10)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=10)
    ax.set_xlim(z.min(), z.max())
    
    if pretty:
        pretty_fig (fig, ax)
    
    # # I used to use this to make the plot window snap to the ticks. 
    # # Probably won't be used again but was a pain in the ass to write.
    
    # # Get major ticks
    # xticks = ax.get_xticks()
    # yticks = ax.get_yticks()
    
    # # Current limits
    # xlim = ax.get_xlim()
    # ylim = ax.get_ylim()
    
    # # Find nearest ticks within range
    # new_xlim = (max([t for t in xticks if t <= xlim[0]]),
    #             min([t for t in xticks if t >= xlim[1]]))
    # new_ylim = (max([t for t in yticks if t <= ylim[0]]),
    #             min([t for t in yticks if t >= ylim[1]]))
    
    # # Apply snapped limits
    # ax.set_xlim(new_xlim)
    # ax.set_ylim(new_ylim)
    
    # def latex_no_trailing_zeros(x, pos):
    #     s = ('%g' % x)
    #     return r'$\mathdefault{%s}$' % s
    
    # ax.xaxis.set_major_formatter(FuncFormatter(latex_no_trailing_zeros))
    # ax.yaxis.set_major_formatter(FuncFormatter(latex_no_trailing_zeros))
    
    return fig, ax

def plot_transparency_mapping_pyvista (mesh : Mesh, arr, cblabel=None, cmap='hot', savelabel='clouds.png', prune=False):
    
    """ REQUIRES pyvista PACKAGE! 
    
    Plots 3D mesh data using transparency mapping. That is, high values = thick cloud, low values = transparent.
    
    Returns an interactive map. Press F2 to screenshot
    
    Parameters
    ----------
    mesh : Mesh
        The mesh in which to solve.
    arr : np.float64 [:]
        Some quantity defined over spatial d.o.f.
    cblabel : str, optional
        Label of the colormap. By default, no label is used.
    cmap : str, optional
        Colormap to use, according to matplotlib
    savelabel : str, optional
        Name of the screenshot file, when F2 is pressed.
    prune : boolean, optional
        Whether or not to take negative values to zero. Further notes on this
        option are in the source code. By default, this is false.
    """
    
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import pyvista as pv
    from Terpdose import Mesh
    
    matplotlib.use('Agg')
    
    # Enable LaTeX
    plt.rcParams.update ({
        'text.usetex' : True,
        'font.family' : 'serif',
        'font.serif'  : ['Computer Modern Roman']
    })
    
    if prune:
        arr[arr < 0] = 0.0
    
    rg = mesh.nodes ()[mesh.connectivity ()]
    o  = mesh.offset ()
    NE = mesh.num_elements
    
    #  ====================================================
    #    Set up the Lionbolt mesh in pyvista's mesh class  
    #  ====================================================
    # MUST BE GENERALIZED TO MORE THAN TETRAHEDRA
    
    cells     = []
    celltypes = []
    
    for e in range(NE):
        ips = o[e]
        ipe = o[e + 1]
        
        nodes = np.arange (ips, ipe)
        
        cells.extend ([4, *nodes])
        celltypes.append (pv.CellType.TETRA)
    
    cells     = np.asarray (cells, dtype=np.int64)
    celltypes = np.asarray (celltypes)
    
    mesh = pv.UnstructuredGrid (cells, celltypes, rg)
    
    # Intensity associated with each discontinuous node
    mesh.point_data[cblabel] = np.asarray(arr)
    
    #  ========
    #    Plot  
    #  ========
    
    p = pv.Plotter(window_size=[1920, 1080])
    
    n = 6
    x = np.linspace(0, 1, n)
    
    gamma = 0.9
    opacity = x**gamma
    
    p.add_volume(
        mesh,
        scalars=cblabel,
        cmap=cmap,
        opacity=opacity,
    )
    
    p.add_axes()
    
    def save_screenshot():
        p.screenshot(savelabel)# , window_size=(1920, 1080)) # , window_size=[1920, 1080])
        print(f"Terpdose.plot_transparency_mapping_pyvista --- Saved screenshot as '{savelabel}'")
    
    p.add_key_event('F2', save_screenshot)
    
    print(f"Terpdose.plot_transparency_mapping_pyvista --- Press 'F2' to save a screenshot as '{savelabel}'")
    print( '                                               (It may take a while depending on the size of your data)')
    
    p.show()
    
    return