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

    # -------------------------------------------------
    # Custom parameters
    # -------------------------------------------------

    #  All examples use u(1) = revs_per_second; u(2) = v_speed
    example = 0               

    # Tfinal
    Tfinal  = 1.0


    # 0 : non-smooth; 
    # 1 : q = 1  (for compressible velocity fields)
    # 2 : smooth (for computing errors)

    init_choice = 0

    # geometry of the cylinder
    R = 1.0
    H = 6.283185307179586

    xc0 = 0.75
    yc0 = 0.25
    r0 = 0.5

    # Speed
    revs_per_s = 1.0
    v_speed = 1.0

    # --------------------- AMR parameters ---------------------------
    # -1 : Always refine
    #  0 : Usual refinement based on threshold
    #  1 : xc < 0.5
    #  2 : yc < 0.5
    refine_pattern = 2

    refine_threshold = 0.05     # Refine everywhere

    dt_initial = 2e-2
    use_fixed_dt = True

    outstyle = 3
    nout = 5
    nstep = 1

    ratioxy = 2
    ratiok = 1

    maxlevel = 2

    # Course grid size
    mx = 16
    my = 16

    exact_metric = 0

    # 0 = no qad
    # 1 = original qad
    # 2 = original (fixed to include call to rpn2qad)
    # 3 = new qad (should be equivalent to 2)
    qad_mode = 0

    fluctuation_sync = 1

    maux = 11
    use_fwaves = True


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('example',          example,        'example')
    probdata.add_param('init_choice',      init_choice,    'init_choice')
    probdata.add_param('refine_pattern',   refine_pattern, 'refine_pattern')
    probdata.add_param('R',                R,              'R')
    probdata.add_param('H',                H,              'H')
    probdata.add_param('xc0',              xc0,            'xc0')
    probdata.add_param('yc0',              yc0,            'yc0')
    probdata.add_param('r0',               r0,             'r0')
    probdata.add_param('revs_per_s',       revs_per_s,     'revs_per_s')
    probdata.add_param('v_speed',          v_speed,        'v_speed')
    probdata.add_param('fluctuation_sync', fluctuation_sync, 'fluctuation_sync')
    probdata.add_param('exact_metric',     exact_metric,    'exact_metric')
  
    probdata.add_param('mx',             mx,             'mx')
    probdata.add_param('mi',             1,              'mi')
    probdata.add_param('mj',             1,              'mj')
    probdata.add_param('maxlevel',       maxlevel,       'maxlevel')
    probdata.add_param('reffactor',      ratioxy,        'reffactor')
    probdata.add_param('qad_mode',       qad_mode,       'qad_mode')

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

    clawdata.output_style = outstyle

    clawdata.dt_variable = not use_fixed_dt
    clawdata.dt_initial = dt_initial

    if clawdata.output_style==1:
        clawdata.num_output_times = nout
        clawdata.tfinal = Tfinal

    elif clawdata.output_style == 2:
        clawdata.output_times =  [0., 0.5, 1.0]

    elif clawdata.output_style == 3:
        clawdata.total_steps = nout
        clawdata.output_step_interval = nstep

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
    clawdata.steps_max = 100000

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 2

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
    clawdata.limiter = ['minmod']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    clawdata.source_split = 0


    # --------------------
    # Boundary conditions:
    # --------------------

    clawdata.num_ghost = 2

    clawdata.bc_lower[0] = 'periodic'   # at xlower
    clawdata.bc_upper[0] = 'periodic'   # at xupper

    clawdata.bc_lower[1] = 'periodic'   # at ylower
    clawdata.bc_upper[1] = 'periodic'   # at yupper


    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    amrdata.amr_levels_max = maxlevel

    amrdata.refinement_ratios_x = [ratioxy]*maxlevel
    amrdata.refinement_ratios_y = [ratioxy]*maxlevel
    amrdata.refinement_ratios_t = [ratiok]*maxlevel

    # If ratio_t == 1
    if ratiok == 1:
        refine_factor = 1
        for i in range(1,maxlevel):
            refine_factor *= amrdata.refinement_ratios_x[i]

        clawdata.dt_initial = dt_initial/refine_factor
        clawdata.total_steps = nout*refine_factor
        clawdata.output_step_interval = nstep*refine_factor

    print("clawdata.dt_initial = {:.6f}".format(clawdata.dt_initial))
    print("clawdata.total_steps = {:d}".format(clawdata.total_steps))
    print("clawdata.output_step_interval = {:d}".format(clawdata.output_step_interval))

    # Refinement threshold
    amrdata.flag2refine_tol = refine_threshold  # tolerance used in this routine


    # ------------------------------------------------------
    # Misc AMR parameters
    # ------------------------------------------------------
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.000000e+00  # Richardson tolerance

    amrdata.flag2refine = True      # use this?

    amrdata.regrid_interval = 100
    amrdata.regrid_buffer_width  = 0
    amrdata.clustering_cutoff = 1 # 0.800000
    amrdata.verbosity_regrid = 0

    # ----------------------------------------------------------------
    # For mapped grid problem
    # 1        capacity
    # 2-5      Velocities projected onto left/right/top/bottom edges
    # 6-9      Edge lengths (normalized by dx or dy)
    # 
    # All values are "center" values, since each cell stores 
    # data for all four faces.  This means that duplicate information
    # is stored, but we have to do it this way for conservation fix.
    # ----------------------------------------------------------------

    amrdata.aux_type = ['capacity'] + ['center']*10

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
