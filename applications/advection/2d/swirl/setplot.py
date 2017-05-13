
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'forestclaw'
    
    # ------------------------------------------------------------
    # Figure for pcolor plot
    # ------------------------------------------------------------

    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[0]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = 0
    plotitem.imshow_cmap = colormaps.yellow_red_blue
    plotitem.imshow_cmin = 0.0
    plotitem.imshow_cmax = 1.0
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [True, True, False]
    plotitem.amr_patchedges_show = [True, True]
    plotitem.show = True       # show on plot?
    
    # ------------------------------------------------------------
    # Figure for contour plot
    # ------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name='contour', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[0]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 0
    plotitem.contour_nlevels = 20
    plotitem.contour_min = 0.01
    plotitem.contour_max = 0.99
    plotitem.amr_contour_colors = ['b','k','r']
    plotitem.show = True       # show on plot?

    # ------------------------------------------------------------
    # Figure for tikz plots (use KML plotting for now)
    # ------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name='swirl (tikz)', figno=2)
    plotfigure.use_for_kml = True
    plotfigure.kml_xlimits = [0,1]
    plotfigure.kml_ylimits = [0,1]

    plotfigure.kml_figsize = [4,4]  #[mx*2**p*mi,mx*2**p*mj]
    plotfigure.kml_dpi = 512

    # Color axis : transparency below 0.1*(cmax-cmin)
    cmin = 0
    cmax = 1
    cmap = colormaps.yellow_red_blue  # transparent --> light blue --> dark blue

    # Water
    plotaxes = plotfigure.new_plotaxes('tikz')
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0   # Plot height field h.    
    plotitem.pcolor_cmap = cmap
    plotitem.pcolor_cmin = cmin
    plotitem.pcolor_cmax = cmax

    def kml_colorbar(filename):
        geoplot.kml_build_colorbar(filename,cmap,cmin,cmax)

    plotfigure.kml_colorbar = kml_colorbar


    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q'

    # Plot q as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'

    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = [2]            # list of figures to print
    plotdata.html = False                    # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.html_movie = False      # new style, or "4.x" for old style
    plotdata.latex = False                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    plotdata.kml = True

    return plotdata

    
