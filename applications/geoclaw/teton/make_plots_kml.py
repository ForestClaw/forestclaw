
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from clawpack.geoclaw import topotools

#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """


    from clawpack.visclaw import colormaps, geoplot
    from numpy import linspace

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = "forestclaw"

    plotdata.verbose = False


    #-----------------------------------------
    # Some global kml flags
    #-----------------------------------------
    plotdata.kml_name = "Teton Dam"
    plotdata.kml_starttime = [1976,6,5,17,55,0]  # Date/time of event in UTC [None]
    plotdata.kml_tz_offset = 6    # Time zone offset (in hours) of event. [None]

    plotdata.kml_index_fname = "TetonDam"  # name for .kmz and .kml files ["_GoogleEarth"]

    # Set to a URL where KMZ file will be published.
    # plotdata.kml_publish = 'http://math.boisestate.edu/~calhoun/visclaw/GoogleEarth/kmz'

    # Add [file_name,visibility]
    plotdata.kml_user_files.append(['teton_dam_validate.kml',True])

    # Cells used in setrun.py (original)
    num_cells = [54,19]
    lower = [-112.36171859324912, 43.591904932832371]
    upper = [-111.25911793671588, 43.977907507732617]
    #xll = [-111.64, 43.913661]
    #xur = [-111.60, 43.92]

    # Lower left   ( -112.34626736,  43.18013542)
    # Upper right  ( -111.26428819,  43.95986458)

    num_cells = [54,54]
    lower = [-112.34626736,  43.18013542]
    upper = [-111.26428819,  43.95986458]

    #xll = [-111.64, 43.913661]
    #xur = [-111.60, 43.92]

    #-----------------------------------------------------------
    # Figure for KML files (large view)
    # This under-resolves the finest level.
    #----------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Teton Dam',figno=1)
    plotfigure.show = True

    plotfigure.use_for_kml = True
    plotfigure.kml_use_for_initial_view = True

    # Latlong box used for GoogleEarth

    plotfigure.kml_xlimits = [lower[0], upper[0]]
    plotfigure.kml_ylimits = [lower[1], upper[1]]

    # Use computational coordinates for plotting
    plotfigure.kml_use_figure_limits = True
    plotfigure.kml_tile_images = False    # Tile images for faster loading.  Requires GDAL [False]


    # --------------------------------------------------
    # Resolution (should be consistent with data)
    # Refinement levels : [2,4,4,4]; max level = 5; num_cells = [55,24]
    # Aim for 1 pixel per finest level grid cell.
    # Choose a figure size that matches the coarse grid resolution.
    # Set the dpi so that figsize*dpi = finest level effective resolution.

    # If amr refinement ratios set to [0,6]; max_level = 6
    # figsize*dpi = [2,1]*16*2**6 = [2048,1024]
    mx = 16
    mi = 1
    mj = 1
    minlevel = 0
    maxlevel = 7
    p = 1
    plotfigure.kml_figsize = [32,32]  #[mx*2**p*mi,mx*2**p*mj]
    plotfigure.kml_dpi = 64

    # --------------------------------------------------

    # Color axis : transparency below 0.1*(cmax-cmin)
    cmin = 0
    cmax = 5
    cmap = geoplot.googleearth_flooding  # transparent --> light blue --> dark blue

    # Water
    plotaxes = plotfigure.new_plotaxes('kml')
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.depth   # Plot height field h.
    plotitem.pcolor_cmap = geoplot.googleearth_flooding
    plotitem.pcolor_cmin = cmin
    plotitem.pcolor_cmax = cmax

    def kml_colorbar(filename):
        geoplot.kml_build_colorbar(filename,cmap,cmin,cmax)

    plotfigure.kml_colorbar = kml_colorbar


    #-----------------------------------------------------------
    # Figure for KML files (zoomed view on region)
    #----------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Power Plant (zoom)',figno=2)
    plotfigure.show = True

    plotfigure.use_for_kml = True
    plotfigure.kml_use_for_initial_view = False

    # Latlong box used for GoogleEarth

    import tools
    #region_lower, region_upper, figsize = tools.region_coords(xll,xur,
    #                                                          num_cells,
    #                                                          lower,
    #                                                          upper)

    # # Get degrees per finest level :
    # dx_coarse = (upper[0]-lower[0])/num_cells[0]
    # dy_coarse = (upper[1]-lower[1])/num_cells[1]
    #
    # # Zoom region
    # mx_xlow = np.floor((xll[0] - lower[0])/dx_coarse).astype(int)
    # my_ylow = np.floor((xll[1] - lower[1])/dy_coarse).astype(int)
    # mx_zoom = np.ceil((xur[0] - xll[0])/dx_coarse).astype(int)
    # my_zoom = np.ceil((xur[1] - xll[1])/dy_coarse).astype(int)
    # xlow = lower[0] + mx_xlow*dx_coarse
    # ylow = lower[1] + my_ylow*dy_coarse

    #plotfigure.kml_xlimits = [region_lower[0],region_upper[0]]
    #plotfigure.kml_ylimits = [region_lower[1], region_upper[1]]

    # Use computational coordinates for plotting
    plotfigure.kml_use_figure_limits = True

    # --------------------------------------------------
    # plotfigure.kml_figsize = figsize*8
    # plotfigure.kml_dpi = 4*4

    # --------------------------------------------------

    plotfigure.kml_tile_images = False    # Tile images for faster loading.  Requires GDAL [False]

    # Color axis : transparency below 0.1*(cmax-cmin)
    # cmin = 0
    # cmax = 5
    # cmap = geoplot.googleearth_flooding  # transparent --> light blue --> dark blue

    # Water
    # plotaxes = plotfigure.new_plotaxes('kml')
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.depth   # Plot height field h.
    # plotitem.pcolor_cmap = geoplot.googleearth_flooding
    # plotitem.pcolor_cmin = cmin
    # plotitem.pcolor_cmax = cmax

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Flood height', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3  # eta?
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = True

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo

    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'

    def afterframe(current_data):
        from pylab import plot, legend, xticks, floor, axis, xlabel,title
        t = current_data.t
        gaugeno = current_data.gaugeno
        if gaugeno == 1:
            title('Wilford')
        elif gaugeno == 2:
            title('Teton City')

        # plot(t, 0*t, 'k')
        n = int(floor(t.max()/3600.) + 2)
        xticks([3600*i for i in range(n)], ['%i' % i for i in range(n)])
        xlabel('time (hours)')

    plotaxes.afteraxes = afterframe


    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.parallel = False
    plotdata.print_format = 'png'           # file format
    plotdata.print_framenos = 'all'         # list of frames to print
    plotdata.print_gaugenos = 'all'         # list of gauges to print
    plotdata.print_fignos = [1,300]         # list of figures to print

    plotdata.printfigs = True              # print figures
    plotdata.overwrite = True

    plotdata.html = False                     # create html files of plots?
    plotdata.html_movie = False                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index

    plotdata.latex = False                    # create latex file of plots?
    #plotdata.latex_figsperline = 2           # layout of plots
    #plotdata.latex_framesperline = 1         # layout of plots
    #plotdata.latex_makepdf = False           # also run pdflatex?

    plotdata.kml = True

    return plotdata

if __name__=="__main__":
    from clawpack.visclaw.plotclaw import plotclaw
    plotclaw(outdir='.',setplot=setplot,plotdir='_plots',format='forestclaw')    
