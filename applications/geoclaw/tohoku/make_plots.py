
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

# import numpy as np
from numpy import *
import matplotlib.pyplot as plt

#-------------------
def setplot(plotdata):
#--------------------


    """
    specify what is to be plotted at each frame
    Input: plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Outptu: a modified version of plotdata.

    """

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=False)

    def fixup(current_data):
        import pylab
        addgauges(current_data)
        t = current_data.t
        t = t / 3600.  # hours
        pylab.title('Surface at %4.2f hours' % t, fontsize=20)
        #pylab.xticks(fontsize=15)
        #pylab.yticks(fontsize=15)


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
    # Some global kml flags
    #-----------------------------------------
    plotdata.kml_index_fname = "Tohoku_2011"     # Name for .kmz and .kml files.
    plotdata.kml_name = "Tohoku 2011"

    plotdata.kml_starttime = [2011,3,11,5,46,0]  # [Y,M,D,H,M,S] (UTC)
    plotdata.kml_tz_offset = -9     # offset to UTC

    plotdata.kml_publish = "http://math.boisestate.edu/~calhoun/visclaw/GoogleEarth/kmz"


    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Bathymetry', figno=0)
    plotfigure.show = True   # Don't show this file in the html version

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    plotaxes.afteraxes = fixup
    plotaxes.xlimits = [132.0, 210.0]
    plotaxes.ylimits = [9.0, 53.0]

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.show = True
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.2
    plotitem.pcolor_cmax = 0.2
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.show = True
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [1,1,0]
    plotitem.patchedges_show = 1
    plotaxes.xlimits = [202., 206.]
    plotaxes.ylimits = [19., 21.]

    # add contour lines of bathy if desired:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.show = True
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = linspace(-2000,0,5)
    plotitem.amr_contour_colors = ['w']  # color on each level
    plotitem.kwargs = {'linestyles':'solid','linewidths':1}
    plotitem.amr_contour_show = [1,0,0]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0



    #-------------------------------------------------------------------
    # Figure for KML files
    #--------------------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Sea Surface',figno=4)
    plotfigure.show = True   # Don't show this file in the html version

    plotfigure.kml_use_for_initial_view= True
    plotfigure.use_for_kml = True

    # Resolution : prod(rr_factors)*22
    rcl = 2    # rcl*figsize = numcells
    mx = 32
    rr_factors = array([1,4,8,4])
    plotfigure.kml_dpi = rr_factors.prod()
    plotfigure.kml_figsize = mx*array([2,1])

    plotfigure.kml_tile_images = False

    # LatLong box used for plotting
    plotfigure.kml_xlimits = [132.0, 210.0]
    plotfigure.kml_ylimits = [9.0, 53.0]

    cmin = -0.2
    cmax = 0.2
    cmap = geoplot.googleearth_transparent

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = cmap
    plotitem.pcolor_cmin = cmin
    plotitem.pcolor_cmax = cmax

    def kml_colorbar(filename):
        geoplot.kml_build_colorbar(filename,
                                   cmap,
                                   cmin,cmax)

    plotfigure.kml_colorbar = kml_colorbar


    #-----------------------------------------
    # Gauge : figure 300
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.axescmd = 'subplot(2,1,1)'
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'
    plotitem.kwargs = {'linewidth':2}

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = False

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo

    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'

    def add_zeroline(current_data):
        from pylab import plot, legend, xticks, floor, xlim,ylim
        t = current_data.t
        plot(t, 0*t, 'k')
        #n = int(floor(t.max()/1800.)) + 2
        #xticks([1800*i for i in range(n)],[str(0.5*i) for i in range(n)])
        #xlim(25000,t.max())
        #ylim(-0.5,0.5)
        print("+++ gaugeno = ",current_data.gaugeno)

    def add_legend_eta(current_data):
        from pylab import legend
        legend(('Surface'),loc='lower left')
        add_zeroline(current_data)

    plotaxes.ylimits = [-1.5, 1.5]
    plotaxes.afteraxes = add_legend_eta


    #-----------------------------------------
    # Gauge : figure 301 (speed, uvel, vvel)
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Velocities', figno=301, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    #plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.title = 'Velocities'
    plotaxes.afteraxes = add_zeroline

    # Plot velocity as red curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = True
    def speed(current_data):
        from numpy import where, sqrt
        h = current_data.q[0,:]
        h = where(h>0.01, h, 1.e6)
        u = 100. * current_data.q[1,:] / h
        v = 100. * current_data.q[2,:] / h
        s = sqrt(u**2 + v**2)
        return s
    plotitem.plot_var = speed
    plotitem.plotstyle = 'k-'

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def uvel(current_data):
        from numpy import where, sqrt
        h = current_data.q[0,:]
        h = where(h>0.01, h, 1.e6)
        u = 100. * current_data.q[1,:] / h
        return u
    plotitem.plot_var = uvel
    plotitem.plotstyle = 'r-'
    plotitem.kwargs = {'linewidth':2}

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    def vvel(current_data):
        from numpy import where, sqrt
        h = current_data.q[0,:]
        h = where(h>0.01, h, 1.e6)
        v = 100. * current_data.q[2,:] / h
        return v
    plotitem.plot_var = vvel
    plotitem.plotstyle = 'g-'
    plotitem.kwargs = {'linewidth':2}

    def add_legend_vel(current_data):
        from pylab import legend
        #legend(['u','v'],loc='upper left')
        add_zeroline(current_data)
        legend(['Speed','u','v'],loc='upper left')

    plotaxes.ylimits = [-50,50]
    plotaxes.afteraxes = add_legend_vel


    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.parallel = False
    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = range(0,101)          # list of frames to print
    plotdata.print_gaugenos = 'all'            # list of gauges to print
    plotdata.print_fignos = [4,300]            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_movie = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = False                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    plotdata.kml = True
    #plotfigure.kml_url = 'http://math.boisestate.edu/~calhoun/visclaw/GoogleEarth/tohoku'

    return plotdata

if __name__=="__main__":
    from clawpack.visclaw.plotclaw import plotclaw
    plotclaw(outdir='.',setplot=setplot,plotdir='_plots',format='forestclaw')    
