"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

from __future__ import absolute_import
import os
import numpy as np

#------------------------------
def setrun(claw_pkg='amrclaw'):
#------------------------------

    from clawpack.clawutil import data


    assert claw_pkg.lower() == 'amrclaw',  "Expected claw_pkg = 'amrclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    # ---------------------------
    # Physical problem parameters
    # ---------------------------

    # 0 = rigid body rotation; 
    # 1 = horizontal flow
    # 2 = sine patch
    # 3 = horizontal flow with variable speed
    example = 3          

    refine_pattern = 1   # 0 = constant theta;  1 = constant_r

    rps   = -1                          # units of theta/second (example=0)
    cart_speed = 1.092505803290319      # Horizontal speed (example=1)
    amplitude = 0.05                    # Amplitude for sine path
    freq = 1                          # Frequency for sine path

    beta = 0.4

    # 0 : top position (12 o'clock)
    # 1 : left position (9 o'clock)
    # 2 : bottom position (6 o'clock)
    # 3 : right position (3 o'clock)
    rotate_position = 0

    flow_direction = 1

    cart_speed = flow_direction*cart_speed

    if rotate_position == 0:
        vcart = [cart_speed,0]
    elif rotate_position == 1:
        vcart = [0,cart_speed]
    elif rotate_position == 2:
        vcart = [-cart_speed, 0]
    else:
        vcart = [0,-cart_speed]    

    theta = np.array([0.125,0.375])  + 0.25*rotate_position
 
    t0 = np.pi/2*(1 + rotate_position + flow_direction*(1/8))
    ravg = (1 + beta)/2
    initial_location = ravg*np.array([np.cos(t0), np.sin(t0)])

    # ---------------
    # Grid parameters
    # ---------------
    grid_mx = 32    # Size of ForestClaw grids
    mi = 4          # Number of ForestClaw blocks
    mj = 2     
    mx = mi*grid_mx
    my = mj*grid_mx

    # -------------
    # Time stepping
    # -------------
    dt_initial = 2.5e-3        # Stable for level 1
    nout = 100                 # 400 steps => T=2
    nsteps = 10

    # ------------------
    # AMRClaw parameters
    # ------------------
    regrid_interval = 10000    # Don't regrid

    maxlevel = 2
    ratioxy = 2
    ratiok = 1

    limiter = 'minmod'    # 'none', 'minmod', 'superbee', 'vanleer', 'mc'

    # 0 = no qad
    # 1 = original qad
    # 2 = original (fixed to include call to rpn2qad)
    # 3 = new qad (should be equivalent to 2)
    qad_mode = 2

    maux = 15
    use_fwaves = True

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

    # example 0 : Rigid body rotation (possibly using a streamfunction)
    # Make vertical speed small so we leave grid
    probdata.add_param('example',                example,           'example')
    probdata.add_param('revolutions per second', rps,               'rps')
    probdata.add_param('vcart[0]',               vcart[0],          'vcart[0]')
    probdata.add_param('vcart[1]',               vcart[1],          'vcart[1]')
    probdata.add_param('amplitude',              amplitude,         'amplitude')
    probdata.add_param('freq',                   freq,              'freq')
    probdata.add_param('initial radius',         0.05,              'init_radius')
    probdata.add_param('initial_location[0]',    initial_location[0], 'initial_location[0]')
    probdata.add_param('initial_location[1]',    initial_location[1], 'initial_location[1]')
    probdata.add_param('beta',                   beta,               'beta')
    probdata.add_param('theta1',                 theta[0],          'theta(1)')
    probdata.add_param('theta2',                 theta[1],          'theta(2)')

    probdata.add_param('grid_mx',        grid_mx,        'grid_mx')
    probdata.add_param('mi',             mi,             'mi')
    probdata.add_param('mj',             mj,             'mj')
    probdata.add_param('maxlevel',       maxlevel,       'maxlevel')
    probdata.add_param('reffactor',      ratioxy,        'reffactor')
    probdata.add_param('refine_pattern', refine_pattern, 'refine_pattern')
    probdata.add_param('rotate_position', rotate_position, 'rotate_position')
    probdata.add_param('qad_mode',        qad_mode,       'qad_mode')

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amrclaw.data for AMR)
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    clawdata.num_dim = num_dim

    clawdata.lower[0] = 0          # xlower
    clawdata.upper[0] = 1          # xupper
    clawdata.lower[1] = 0          # ylower
    clawdata.upper[1] = 1          # yupper

    clawdata.num_cells[0] = mx      # mx
    clawdata.num_cells[1] = my      # my

    clawdata.num_eqn = 1
    clawdata.num_aux = maux
    clawdata.capa_index = 1


    # ----------------------------------------------------------
    # Time stepping
    # ----------------------------------------------------------

    clawdata.output_style = 3

    clawdata.dt_variable = False
    clawdata.dt_initial = dt_initial

    if clawdata.output_style==1:
        clawdata.num_output_times = 16
        clawdata.tfinal = 4.0

    elif clawdata.output_style == 2:
        clawdata.output_times =  [0., 0.5, 1.0]

    elif clawdata.output_style == 3:
        clawdata.total_steps = nout
        clawdata.output_step_interval = nsteps

    clawdata.output_format = 'ascii'       # 'ascii', 'binary', 'netcdf'

    # ---------------------------
    # Misc time stepping and I/O
    # ---------------------------
    clawdata.cfl_desired = 0.900000
    clawdata.cfl_max = 1.000000

    clawdata.output_t0 = True  # output at initial (or restart) time?
    clawdata.t0 = 0.000000

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data


    clawdata.output_q_components = 'all'    # only 'all'
    clawdata.output_aux_components = 'none'  # 'all' or 'none'
    clawdata.output_aux_onlyonce = False    # output aux arrays only at t0?


    clawdata.dt_max = 1.000000e+99
    clawdata.steps_max = 1000

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = maxlevel

    # ----------------------------------------------------
    # Clawpack parameters
    # -----------------------------------------------------

    clawdata.order = 2
    clawdata.dimensional_split = 'unsplit'
    clawdata.transverse_waves = 2

    clawdata.num_waves = 1

    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = [limiter]

    clawdata.use_fwaves = use_fwaves    # True ==> use f-wave version of algorithms

    clawdata.source_split = 0


    # --------------------
    # Boundary conditions:
    # --------------------

    clawdata.num_ghost = 2

    clawdata.bc_lower[0] = 'extrap'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper

    clawdata.bc_lower[1] = 'extrap'   # at ylower
    clawdata.bc_upper[1] = 'extrap'   # at yupper


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    amrdata.amr_levels_max = maxlevel

    amrdata.refinement_ratios_x = [ratioxy]*maxlevel
    amrdata.refinement_ratios_y = [ratioxy]*maxlevel
    amrdata.refinement_ratios_t = [ratiok]*maxlevel

    # If we are taking a global time step (stable for maxlevel grids), we 
    # probably want to increase number of steps taken to hit same
    # time 'tfinal'.  
    if ratiok == 1:
        refine_factor = 1
        for i in range(1,maxlevel):
            refine_factor *= amrdata.refinement_ratios_x[i]

        # Decrease time step
        clawdata.dt_initial = dt_initial/refine_factor

        # Increase number of steps taken.
        clawdata.total_steps = nout*refine_factor
        clawdata.output_step_interval = nsteps*refine_factor

    # Refinement threshold
    amrdata.flag2refine_tol = -1  # tolerance used in this routine


    # ------------------------------------------------------
    # Misc AMR parameters
    # ------------------------------------------------------
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.000000e+00  # Richardson tolerance

    amrdata.flag2refine = True      # use this?

    amrdata.regrid_interval = regrid_interval
    amrdata.regrid_buffer_width  = 0
    amrdata.clustering_cutoff = 0.800000
    amrdata.verbosity_regrid = 0


    # ----------------------------------------------------------------

    # Conservative form (cell-centered velocities)
    # 1      capacity
    # 2-3    Cell-centered velocities projected 
    # 4-5    normal at x face
    # 6-7    normal at y face
    # 8-9    edgelengths at x/y faces
    # ----------------------------------------------------------------

    if qad_mode in [0,1]:
        amrdata.aux_type = ['capacity'] + ['center']*2 + ['xleft']*4 + \
        ['yleft']*4 + ['xleft']*2 + ['yleft']*2
    else:
        # Each cell has data for all four faces
        amrdata.aux_type = ['capacity'] + ['center']*14

    #  ----- For developers -----
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False       # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting

    return rundata

    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
