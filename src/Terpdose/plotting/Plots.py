from __future__ import annotations

# POSSIBLE PLAN --- Make these functions accept mesh, array, geo,
# then they figure out how to actually plot.
# Or just make a 1D plot 2D plot etc.
# Could still have wrapper functions that 

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

def plot_1D (mesh, geo : Line, array, title=None, xlabel=r'$d$ (cm)', ylabel=None, color='blue', pretty=True):
    
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
        at the very end by sending fig and ax to the function pretty_fig (fig, ax)
    
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
    plt.rcParams.update({
        'text.usetex' : True,
        'font.family' : 'serif',
        'font.serif'  : ['Computer Modern Roman']
    })
    
    u = np.array([np.dot(geo.xyz[i,:] - geo.x0, geo.ehat) for i in range(geo.n)])
    y = GeoInterpolation (mesh, geo, array)
    
    fig, ax = plt.subplots(figsize=(7.1, 4))
    
    ax.plot (u, y, '-o', markersize=0.5, linewidth=0.75, color=color)
    
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
    plt.rcParams.update({
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

def plot_slab (mesh : Mesh, array, FMR=False, title=None, xlabel=r'Depth (cm)', ylabel=None, color='blue', pretty=True):
    
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
        at the very end by sending fig and ax to the function pretty_fig (fig, ax)
    
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
    plt.rcParams.update({
        'text.usetex' : True,
        'font.family' : 'serif',
        'font.serif'  : ['Computer Modern Roman']
    })
    
    z = mesh.extended_nodes()[:,0]
    
    if FMR:
        z = z / max(z)
    
    fig, ax = plt.subplots(figsize=(7.1, 4))
    
    ax.plot (z, array, '-o', markersize=0.5, linewidth=0.75, color=color)
    
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
