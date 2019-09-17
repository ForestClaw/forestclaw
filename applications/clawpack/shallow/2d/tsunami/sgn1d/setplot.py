
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

from numpy import *
#from matplotlib.pylab import *

# import pyclaw.data

#--------------------------
def setplot(plotdata):
#--------------------------

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    # -------------------------------------------------------------------------
    # Plot figures
    # -------------------------------------------------------------------------

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Height', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-100,500]
    plotaxes.ylimits = [0.975,1.15]
    plotaxes.title = 'Height'

    # Set up for item on these axes:
    # plotaxes.afteraxes = plot_exact_1

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    # plotitem.show = True       # show on plot?


#     # Plot linear solution
#     plotitem = plotaxes.new_plotitem(name='line2',plot_type='1d_plot')
#     plotitem.outdir = '../../../ex_linear_swe/linear_swe/_output'
#     plotitem.plot_var = 0
#     plotitem.plotstyle = '-'
#     plotitem.color = 'r'
#     plotitem.show = True       # show on plot?


    plotfigure = plotdata.new_plotfigure(name='Velocity', figno=2)

    plotaxes = plotfigure.new_plotaxes()
    # plotaxes.axescmd = 'subplot(2,1,2)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-1, 3]
    plotaxes.title = 'Solution q(2)'

    #plotaxes.afteraxes = plot_exact_2

    # Set up for item on these axes:
    # q = current_data.q
    # h = q[0]
    # u = q[1]/q[0]

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 1
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    # plotitem.show = True       # show on plot?


    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = [1]            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = False                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata
