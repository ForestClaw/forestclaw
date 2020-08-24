"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import numpy as np
from pdb import *

# import topotools
try:
    FCLAW = os.environ['FCLAW']
except:
    raise Exception("*** Must First set FCLAW environment variable")

topodir = os.path.join(FCLAW, 'applications','geoclaw', 'scratch')
if not os.path.isdir(topodir):
    raise Exception("*** Missing topo directory: %s" % topodir)

minlevel = 3
maxlevel = 7



#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2

    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Parameters
    # ---------------

    mx = 32    # Ignored
    my = 32




    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = 132.0          # xlower
    clawdata.upper[0] = 210.0          # xupper
    clawdata.lower[1] = 9.0          # ylower
    clawdata.upper[1] = 53.0          # yupper

    # Number of grid cells (ignored by ForestClaw)
    clawdata.num_cells[0] = mx      # mx
    clawdata.num_cells[1] = my      # my


    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # ---------------
    # Gauges:
    # ---------------
    gauges = rundata.gaugedata.gauges

    # for gauges append lines of the form  [gaugeno, x, y, t1, t2] (25200,inf)
    gauges.append([1123, 203.52825, 20.9021333, 7.0*3600., 1.e9]) #Kahului

    #gauges.append([5680, 203.52333, 20.895, 7.0*3600., 1.e9]) #TG Kahului
    # more accurate coordinates from Yong Wei at PMEL:  (25200,inf)
    gauges.append([5680, 203.530944, 20.895, 7.0*3600., 1.e9]) #TG Kahului


    ## Add gauges for comparison with GeoClaw
    gauges.append([1, 143.214480, 38.011905, 0, 1.e9])
    gauges.append([2, 144.851543, 33.090886, 0, 1.e9]) 
    gauges.append([3, 170.233553, 29.893284, 0, 1.e9])
    gauges.append([4, 196.417438, 20.561113, 0, 1.e9])


    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 26
        clawdata.tfinal = 13*3600.
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        clawdata.output_times =  np.linspace(7,13,4)*3600.

    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = False  # output at initial (or restart) time?

    clawdata.output_format = 'ascii'      # 'ascii' or 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.016

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 50000




    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2

    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'

    # For unsplit method, transverse_waves can be
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif clawdata.checkpt_style == 2:
        # Specify a list of checkpoint times.
        clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5


    # -----------------------------------------------
    # AMR parameters:
    # -----------------------------------------------
    amrdata = rundata.amrdata

    maxlevel = 15

    # List of refinement ratios at each level (length at least amr_level_max-1)
    # 2 degree, 24', 4', 1', 10", 1/3"
    amrdata.amr_levels_max = maxlevel    # Set to 3 for best results

    # Not used in ForestClaw
    amrdata.refinement_ratios_x = [2]*(maxlevel-1)
    amrdata.refinement_ratios_y = [2]*(maxlevel-1)
    amrdata.refinement_ratios_t = [2]*(maxlevel-1)


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?

    amrdata.flag2refine = True
    amrdata.flag2refine_tol = 0.5  # tolerance used in this routine

    amrdata.regrid_interval = 3
    amrdata.regrid_buffer_width  = 2
    amrdata.clustering_cutoff = 0.700000
    amrdata.verbosity_regrid = 0


    # ---------------
    # Regions:
    # ---------------
    regions = rundata.regiondata.regions

    # ------------------------------
    # Regions
    # ------------------------------

    inf = 1e9

    # Not sure why these are included : 
    # regions.append([minlevel, minlevel+3, 0., 5.*3600., 132., 220., 5., 40.])

    # Region describing topo region
    # Data on topo extent extracted from Fujii.txydz
    # 
    # ------------------------------------------------------------------
    # To match Geoclaw levels (N0 = mx = 16)
    # 
    #            GeoClaw              -->           ForestClaw
    #
    #  Levels   RR    RR.cumprod()             Levels   RR   RR.cumprod()     
    #   1       (1)        (1)         -->       0     (1)        (1)
    #   2       (5)        (5)         -->       2     (4)        (4)
    #   3       (6)       (30)         -->       5     (8)       (32)
    #   4       (4)      (120)         -->       7     (4)      (128)
    #   5       (6)      (720)         -->       9     (4)      (512)
    #   6      (30)    (21600)         -->       14   (32)    (16384)
    # ------------------------------------------------------------------

    # Region 0  : Encompasses entire domain (limits refinement 
    # in upper 1/4 portion of domain)  (0, inf)
    regions.append([0, 1, 0., 1e9, 0, 360, -90, 90])

    # Region 1  : Topo map at initial earthquake site;  Refine long enough to 
    # resolve initial disturbance.  (0,3600)
    regions.append([5, 5, 0., 1, 135., 150., 30., 45.])

    # Region 2 : Large region encompassing most of lower portion of domain (0,18000)
    regions.append([0, 5, 0., 5.*3600., 132., 220., 5., 40.])

    # Region 3 :  Large region encompassing Hawaii Islands (18000,28800)
    regions.append([0, 5, 5.0*3600.,  8.0*3600, 180.0, 220.0,  5.0, 40.0])

    # Region 4  : Includes Maui and Molokai (23400.0, inf)
    regions.append([7, 7, 6.5*3600.,  inf,      202.5, 204.0, 20.4, 21.4])

    # Region 5  : Strip including north shore of Maui  (25200, inf)
    regions.append([9, 9, 7.*3600., inf, 203.0, 203.7, 20.88333, 21.])

    # Region 6 : Port at Kailua  (26100.0, inf)
    regions.append([14, 14, 7.25*3600., inf, 203.52,203.537,20.89,20.905])

    # ------------------------------------------------------------------
    # Try to reduce amount of time spent in ghost filling
    # To match Geoclaw levels (mx = 32)
    # 
    #            GeoClaw              -->           ForestClaw
    #
    #  Levels   RR      RR.cumprod()         Levels    RR      RR.cumprod()
    #   1       (1)        (1)         -->       0     (1)        (1)
    #   2       (5)        (5)         -->       2     (4)        (4)
    #   3       (6)       (30)         -->       5     (8)       (32)
    #   4       (4)      (120)         -->       6     (2)       (64)   < ---- only difference
    #   5       (6)      (720)         -->       9     (8)      (512)
    #   6      (30)    (21600)         -->      14    (32)    (16384)
    # ------------------------------------------------------------------

#    # Region 0  : Encompasses entire domain (limits refinement 
#    # in upper 1/4 portion of domain)  (0, inf)
#    regions.append([0, 1, 0., 1.0e9, 0, 360, -90, 90])
#
#    # Region 1  : Topo map at initial earthquake site;  
#    # Extent of the region (including time component) is taken from the dtopo file.
#    regions.append([5, 5, 0, 1, 135., 150., 30., 45.])
#
#    # Region 2 : Large region encompassing most of lower portion of domain (0,18000)
#    regions.append([0, 5, 0, 5*3600., 132., 220., 5., 40.])
#
#    # Region 3 :  Large region encompassing Hawaii Islands (18000,28800)
#    regions.append([0, 5, 5.0*3600.,  8.0*3600, 180.0, 220.0,  5.0, 40.0])
#
#    # Region 4  : Includes Maui and Molokai (23400.0, inf)
#    regions.append([6, 6, 6.5*3600.,  inf,      202.5, 204.0, 20.4, 21.4])
#
#    # Region 5  : Strip including north shore of Maui  (25200, inf)
#    regions.append([9, 9, 7.*3600., inf, 203.0, 203.7, 20.88333, 21.])
#
#    # Region 6 : Port at Kailua  (26100.0, inf)
#    regions.append([14, 14, 7.25*3600., inf, 203.52,203.537,20.89,20.905])


    # -------------------------------------------------------
    # For developers
    #    -- Toggle debugging print statements:
    # -------------------------------------------------------
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = True      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting


    return rundata
    # end of function setrun
    # ----------------------


#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")


    topofile = 'topos/TetonLarge.topo'

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2   # LatLong coordinates
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.035
    geo_data.friction_depth = 500

    # Refinement data
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.016    # Original setting : 0.016
    refinement_data.deep_depth = 200
    refinement_data.max_level_deep = 4

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]

    # == settopo.data values ==
    # Set minlevel=maxlevel=0
    topofiles = rundata.topo_data.topofiles
    topofiles.append([3, 0, 0, 0.0, 1e10, os.path.join(topodir,'etopo1min130E210E0N60N.asc')])
    # topofiles.append([3, 0, 0, 0.0, 1e10, os.path.join(topodir,'hawaii_6s.txt')])
    topofiles.append([3, 0, 0, 0., 1.e10, os.path.join(topodir,'kahului_1s.txt')])

    # == setdtopo.data values ==
    # topo_data = rundata.topo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, minlevel,maxlevel,fname]
    rundata.dtopo_data.dtopofiles = [[1, 0, 0,os.path.join(topodir,'fujii.txydz')]]


    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    # == fixedgrids.data values ==
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    rundata.fixed_grid_data.fixedgrids = []
    fixedgrids = rundata.fixed_grid_data.fixedgrids

    return rundata
    # end of function setgeo
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
