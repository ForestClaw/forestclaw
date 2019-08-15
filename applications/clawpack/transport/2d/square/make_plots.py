
""" 
Plot swirl using Clawpack's visclaw graphics.  This file can be run as ; 

    % python make_plots.py

To learn more about visclaw graphics, see www.clawpack.org
    
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
    import clawpack.forestclaw as pyclaw

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'forestclaw'
    
    example = 1

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
    plotitem.imshow_cmax = 16.0
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [False, False, False]
    plotitem.amr_patchedges_show = [True, True]
    plotitem.show = True       # show on plot?
    
    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = [0]            # list of figures to print
    plotdata.html = True                    # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.html_movie = 'JSAnimation'      # new style, or "4.x" for old style
    plotdata.latex = False                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

if __name__=="__main__":
    from clawpack.visclaw.plotclaw import plotclaw
    plotclaw(outdir='.',setplot=setplot,plotdir='_plots',format='forestclaw')



    
