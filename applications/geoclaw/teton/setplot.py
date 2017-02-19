
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""
"""
Set up the plot figures, axes, andd items to be done for each frame.

this module is imported by the plotting routines and then the
fucntion setplot is called to set the plot parameters.

"""

import numpy as np
import matplotlib.pyplot as plt

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


    #to plot gauge locations on pcolor or contour, use this as
    #after axis function:


import numpy
#a = 1.
#sigma = 0.5
#h0 = 150
#grav = 9.81

#-------------------
def setplot(plotdata):
#--------------------


    """
    specify what is to be plotted at each frame
    Input: plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Outptu: a modified version of plotdata.

    """



    from clawpack.visclaw import colormaps, geoplot

    plotdata.clearfigures()

    def set_drytol(current_data):
        # The drytol parameter is used in masking land and water and
        # affects what color map is used for cells with small water depth h.
        # The cell will be plotted as dry if h < drytol.
        # The best value to use often depends on the application and can
        # be set here (measured in meters):
        current_data.user["drytol"] = 1.e-3

    plotdata.beforeframe = set_drytol

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Teton Dam', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True


    dark_blue = [0.2,0.2,0.7];
    light_blue = [0.7,0.7,1.0];
    flooding_colormap = colormaps.make_colormap({ -1.0:light_blue,
                                                 1.0:dark_blue})
    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface
    plotitem.pcolor_cmap = flooding_colormap
    plotitem.pcolor_cmin = 0
    plotitem.pcolor_cmax = 30
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = False

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 900.0
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 0
    plotaxes.xlimits = [0,48000]
    plotaxes.ylimits = [0,17400]

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = linspace(0,900,40)
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True


    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
                                        gaugenos='5', format_string='ko', add_labels=True)

    plotaxes.afteraxes = addgauges

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
    plotitem.plot_var = 3
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

    plotdata.parallel = True
    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = range(0,36)          # list of frames to print
    plotdata.print_gaugenos = 'all'            # list of gauges to print
    plotdata.print_fignos = [0,300]            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_movie = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = False                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.kml = False

    return plotdata
