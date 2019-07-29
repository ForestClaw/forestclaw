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

    refine_threshold = 15     # Refine everywhere
    dt_initial = 4e-3          # Stable for level 1
    use_fixed_dt = True

    outstyle = 3
    nout = 30
    nstep = 1

    example = 0
    alpha = 0.4
    beta = 0


    maxlevel = 3
    ratioxy = 4
    ratiok = 1


    grid_mx = 32
    mi = 5
    mj = 2
    mx = mi*grid_mx
    my = mj*grid_mx


    # 0 = no qad
    # 1 = original qad
    # 2 = original (fixed to include call to rpn2qad)
    # 3 = new qad (should be equivalent to 2)
    qad_mode = 3

    maux = 9
    use_fwaves = True


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    # c      example_in,mapping_in, ic_in, alpha_in,rps_in)

    probdata.add_param('example',        example,        'example')
    probdata.add_param('alpha',          alpha,          'alpha')
    probdata.add_param('beta',           beta,           'beta')

    probdata.add_param('grid_mx',        grid_mx,        'grid_mx')
    probdata.add_param('mi',             mi,             'mi')
    probdata.add_param('mj',             mj,             'mj')
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
        clawdata.num_output_times = 16
        clawdata.tfinal = 4.0

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
    clawdata.limiter = ['vanleer']

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

    amrdata.regrid_interval = 3
    amrdata.regrid_buffer_width  = 3
    amrdata.clustering_cutoff = 0.800000
    amrdata.verbosity_regrid = 0

    # ----------------------------------------------------------------
    # For torus problem
    # 1        capacity
    # 2-3      Edge velocities
    # 4-6      Cell-centered velocities (x,y,z)
    # 7-8      Edge lengths (x-face, y-face)
    # 9-11     x-face normals
    # 12-14    y-face normals
    # ----------------------------------------------------------------

    if (qad_mode in [0,1]):
        amrdata.aux_type = ['capacity'] + ['xleft','center','yleft','center']*2
    else:
        amrdata.aux_type = ['capacity'] + ['center']*8

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
