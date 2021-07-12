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

    refine_threshold = 0.05     # Refine everywhere
    use_fixed_dt = True

    mx = 32
    dt_initial = 5e-3

    # -1 : Always refine
    #  0 : Usual refinement based on threshold
    #  1 : constant theta
    #  2 : constant r
    refine_pattern = 1

    f = 5

    uniform = True
    if uniform:
        # No refinement;  just increase resolution on each patch
        maxlevel = 1  
        mx *= 2**f    # Coarsest level
        dt_initial /= 2**f
        nout = 5*2**f
        nstep = 5*2**f


    adapt = not uniform
    if adapt:
        # Refinement everywhere
        refine_pattern = -1   # always adapt
        refine_pattern = 0   # usual refinement
        maxlevel = f + 1     # refine everywhere to this level
        nout = 50            # Number of coarse grid steps
        nstep = 50

    outstyle = 3

    # 0 : Rigid body rotation
    # 1 : Vertical flow u = (0,omega)
    # 2 : Horizontal flow (doesn't work for the torus)
    # 3 : Swirl example
    example = 3


    # 0 : non-smooth; 
    # 1 : q = 1  (for compressible velocity fields)
    # 2 : smooth (for computing errors)

    init_choice = 2

    alpha = 0.4
    beta = 0

    # theta_range = [0.125, 0.375]
    theta_range = [0,1]
    phi_range = [0,1]

    init_radius = 0.4

    if example == 0:
        revs_per_s = 1.0
    elif example >= 1:
        revs_per_s = 1.0

    cart_speed = 0.765366864730180 

    ratioxy = 2
    ratiok = 2

    grid_mx = mx
    mi = 5
    mj = 2
    mx = mi*grid_mx
    my = mj*grid_mx


    # 0 = no qad
    # 1 = original qad
    # 2 = original (fixed to include call to rpn2qad)
    # 3 = new qad (should be equivalent to 2)
    qad_mode = 2

    maux = 11
    use_fwaves = True


    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    probdata.add_param('example',        example,        'example')
    probdata.add_param('init_choice',    init_choice,    'init_choice')
    probdata.add_param('refine_pattern', refine_pattern, 'refine_pattern')
    probdata.add_param('alpha',          alpha,          'alpha')
    probdata.add_param('beta',           beta,           'beta')
    probdata.add_param('init_radius',    init_radius,    'init_radius')
    probdata.add_param('revs_per_s',     revs_per_s,     'revs_per_s')
    probdata.add_param('cart_speed',     cart_speed,     'cart_speed')
    probdata.add_param('theta0',         theta_range[0], 'theta[0]')
    probdata.add_param('theta1',         theta_range[1], 'theta[1]')
    probdata.add_param('phi0',           phi_range[0],   'phi[0]')
    probdata.add_param('phi1',           phi_range[1],   'phi[1]')
  
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
    clawdata.verbosity = 1

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

    amrdata.regrid_interval = 1
    amrdata.regrid_buffer_width  = 3
    amrdata.clustering_cutoff = 0.800000
    amrdata.verbosity_regrid = 0

    # ----------------------------------------------------------------
    # For torus problem
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
