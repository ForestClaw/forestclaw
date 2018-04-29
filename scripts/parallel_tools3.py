import os
import sys
import numpy as np
from string import Template
import re
import configparser
from pdb import set_trace


def compile_results(results_dir=None,results_file='results.out',
                    execname=None,duplicate=False):

    if results_dir == None:
        results_dir = os.getcwd()

    # Stuff to grep for
    stats_list_float = [ 'INIT',
                         'ADVANCE$',
                         'GHOSTFILL$',
                         'REGRID$',
                         'GHOSTPATCH_COMM',
                         'ADAPT_COMM',
                         'PARTITION_COMM',
                         'CFL_COMM',
                         'WALLTIME',
                         'GHOSTFILL_COPY',
                         'LOCAL',
                         'COMM',
                         'PARTITION',
                         'GHOSTPATCH_BUILD',
                         'GHOSTFILL_STEP1',
                         'GHOSTFILL_STEP2',
                         'GHOSTFILL_STEP3',
                         'GHOSTFILL_COPY',
                         'GHOSTFILL_INTERP',
                         'GHOSTFILL_AVERAGE',
                         'LOCAL_BOUNDARY_RATIO',
                         'REMOTE_BOUNDARY_RATIO']

    stats_list_int = ['ADVANCE_STEPS_COUNTER$',
                      'GRIDS_PER_PROC$',
                      'LOCAL_BOUNDARY',
                      'REMOTE_BOUNDARY',
                      'INTERIOR']

    # ------------------------------
    # Options : Variable field width
    # ------------------------------
    option_list = {'jobid'  :('jobid', '{jobid[0]:>6s}',   '{:>6d}'),
                   'p'      :('p',     '{p[0]:>8s}',       '{:>8d}'),
                   'mx'     :('mx',    '{mx[0]:>6s}',      '{:>6d}'),
                   'mi'     :('mi',    '{mi[0]:>4s}',      '{:>4d}'),
                   'mj'     :('mj',    '{mj[0]:>4s}',      '{:>4d}'),
                   'min'    :('min',   '{min[0]:>6s}',     '{:>6d}'),
                   'max'    :('max',   '{max[0]:>6s}',     '{:>6d}'),
                   'nout'   :('nout',  '{nout[0]:>8s}',    '{:>8d}'),
                   'tfinal' :('tfinal','{tfinal[0]:>12s}', '{:>12.2e}')}


    option_list_header = "{jobid[1]:>6s}{p[1]:>8s}{mx[1]:>6s}{mi[1]:>4s}{mj[1]:>4s}" + \
                         "{min[1]:>6s}{max[1]:>6s}{nout[1]:>8s}{tfinal[1]:>12s}"
    option_fmt_header = option_list_header.format(**option_list)  # Process names

    option_str_header = option_fmt_header.format(**option_list)
    option_width = 60     # add up all fields widths, i.e. 6+8+6+4+... = 60

    # ---------------------
    # Ints : Field width 12
    # ---------------------
    int_list = ('adv_steps','grids_proc','local_bdry','remote_bdry','interior')
    int_width = 12*len(int_list)
    int_str = ("{:>12s}"*len(int_list)).format(*int_list)

    # -----------------------
    # Floats : Field width 12
    # -----------------------
    float_list = ('init','advance','ghostfill','regrid','patch_comm',
                  'adapt','partition','cfl','walltime',
                  'gf_copy','local','comm','partition','gbuild','step1',
                  'step2','step3', 'copy',
                  'interp','average','l_ratio','r_ratio')
    float_width = 12*len(float_list)
    float_str = ("{:>12s}"*len(float_list)).format(*float_list)


    # Compile data to this file
    line_len = option_width + int_width + float_width
    header_str = option_str_header + int_str + float_str

    # Write out header to file
    resultsfile = open(results_file,'w')
    # resultsfile.write("# " + "-"*line_len)
    # resultsfile.write("\n")

    resultsfile.write(header_str)

    # resultsfile.write("# " + "-"*line_len)
    resultsfile.write("\n")

    # Look for files that look like "torus_00004.o4567".  But if execname is
    # not supplied, be happy with "*_00004.o4567"
    if not execname == None:
        pattern = re.compile("%s_[0-9]{5}\.o[0-9]*" % execname)
    else:
        pattern = re.compile(".*_[0-9]{5}\.o[0-9]*")

    output_files = os.listdir(results_dir)
    for f in output_files:
        if re.match(pattern,f):
            # Extract proc count and jobid from file name, e.g. torus_00016.o45678
            s = f.partition('_')[2].partition('.o')

            pcount = int(s[0])
            jobid = int(s[2])

            run_file = open(f,'r')
            lines = run_file.readlines()

            found_bad_file = False

            for i,l in enumerate(lines):
                if re.search("exceeded",l):
                    found_bad_file = True
                    break

            #if found_bad_file:
            #    print "Skipping file %s (time exceeded)" % f
            #    continue

            # jobid
            resultsfile.write("  ")   # Make up for "# " added as comments above
            s = "{jobid[2]:s}".format(**option_list)
            resultsfile.write(s.format(jobid))

            # proc count
            resultsfile.write("%8d" % pcount)

            # Get values of options

            # mx
            for i,l in enumerate(lines):
                found_mx = False
                if re.search("mx",l):
                    l1 = lines[i].split()
                    mx = int(l1[2])
                    s = "{mx[2]:s}".format(**option_list)
                    resultsfile.write(s.format(mx))
                    found_mx = True
                    break
                if not found_mx:
                    if re.search("nxmax",l):
                        l1 = lines[i].split()
                        mx = int(l1[2])
                        s = "{mx[2]:s}".format(**option_list)  # for ash3d model
                        resultsfile.write(s.format(mx))
                        found_mx = True
                        break


            # mi (look for something distinctive;  count down from there)
            mi = 1
            mj = 1
            for i,l in enumerate(lines):
                if re.search("manifold",l):  # 'mi' not distinctive enough
                    if not re.search("mi",lines[i+1]):
                        print("Brick count mi not found in {:s}; setting mi=1".format(f))
                        mi = 1
                    else:
                        mi_line = lines[i+1].split()
                        mi = int(mi_line[2])
                    if not re.search("mj",lines[i+2]):
                        print("Brick count mj not found in {:s}; setting mj=1".format(f))
                        mj = 1
                    else:
                        mj_line = lines[i+2].split()
                        mj = int(mj_line[2])

                    s = "{mi[2]:s}{mj[2]:s}".format(**option_list)
                    resultsfile.write(s.format(mi,mj))
                    break
                    # print "mi = %d; mj = %d" %(mi,mj)

            # minlevel
            for i,l in enumerate(lines):
                if re.search("minlevel",l):
                    l1 = lines[i].split()
                    if duplicate:
                        minlevel = int(l1[2]) + np.log(float(mi))/np.log(2.0)
                    else:
                        minlevel = int(l1[2])

                    s = "{min[2]:s}".format(**option_list)
                    resultsfile.write(s.format(int(minlevel)))
                    break

            # maxlevel
            for i,l in enumerate(lines):
                if re.search("[^-]maxlevel",l):  # don't match maxlevel-uses-outstyle
                    l1 = lines[i].split()
                    if duplicate:
                        maxlevel = int(l1[2]) + np.log(float(mi))/np.log(2.0)
                    else:
                        maxlevel = int(l1[2])

                    s = "{max[2]:s}".format(**option_list)
                    resultsfile.write(s.format(int(maxlevel)))
                    break

            # nout
            for i,l in enumerate(lines):
                if re.search("nout",l):
                    l1 = lines[i].split()
                    nout = int(l1[2])
                    s = "{nout[2]:s}".format(**option_list)
                    resultsfile.write(s.format(int(l1[2])))
                    break

            # tfinal
            for i,l in enumerate(lines):
                if re.search("tfinal",l):
                    l1 = lines[i].split()
                    tfinal = float(l1[2])
                    s = "{tfinal[2]:s}".format(**option_list)
                    resultsfile.write(s.format(float(l1[2])))
                    break


            # Get counter values
            for k,w in enumerate(stats_list_int):
                found = False
                for i,l in enumerate(lines):
                    if re.search(w,l):
                        a = i
                        found = True
                        break

                if found:
                    l2 = lines[a+2].split()
                    resultsfile.write("%12d" % np.round(float(l2[5])).astype(int))
                else:
                    print("{:s} in file {:s} not found".format(w,f))
                    resultsfile.write("%12f" % (np.nan))

            # Get timer values
            for k,w in enumerate(stats_list_float):
                found = False
                for i,l in enumerate(lines):
                    if re.search(w,l):
                        a = i
                        found = True
                        break

                if found:
                    l2 = lines[a+2].split()
                    resultsfile.write("%12.4e" % float(l2[5]))
                else:
                    print("{:s} in file {:s} not found".format(w,f))
                    resultsfile.write("%12f" % (np.nan))

            resultsfile.write('\n')

            run_file.close()

    resultsfile.close()



def write_ini_files(input_file='create_run.ini',problem='advection'):
    config = configparser.SafeConfigParser(allow_no_value=True)
    config.read(input_file)

    scheduler     = config.get('Defaults','scheduler').partition('#')[0].strip()


    # Example number
    try:
        example    = int(config.get('Defaults','example').partition('#')[0].strip())
    except:
        example = 0


    # Name to use in Juqueen output file
    execname      = config.get('Problem','execname').partition('#')[0].strip()

    # Estimated time per grid
    time_per_grid = float(config.get('Problem','time_per_grid').partition('#')[0].strip())

    # Fit to model
    try:
        coeff_A = float(config.get('Problem','adapt_coeff_A').partition('#')[0].strip())
        coeff_B = float(config.get('Problem','adapt_coeff_B').partition('#')[0].strip())
    except:
        coeff_A = 1.0
        coeff_B = 0


    # Brick dimensions
    try:
        mi0 = int(config.get('Run','mi').partition('#')[0].strip())
        mj0 = int(config.get('Run','mj').partition('#')[0].strip())
    except:
        mi0 = 1
        mj0 = 1


    mx0     = int(config.get('Run', 'mx').partition('#')[0].strip())
    minlevel0  = int(config.get('Run', 'minlevel').partition('#')[0].strip())
    maxlevel0  = int(config.get('Run', 'maxlevel').partition('#')[0].strip())
    proc0   = int(config.get('Run','proc').partition('#')[0].strip())
    try:
        tfinal0 = float(config.get('Run','tfinal').partition('#')[0].strip())
    except:
        pass


    # nbjobs : determines length of processor sequence, e.g. [1,4,16,64,...]
    njobs   = int(config.get('Run','njobs').partition('#')[0].strip())

    if minlevel0 == maxlevel0:
        mode = 'uniform'

    try:
        # Set regrid here in case we want to regrid even in the
        # uniform case
        regrid_interval = int(config.get('Run','regrid_interval').partition('#')[0].strip())
    except:
        regrid_interval = 1




    # Ranks per node
    try:
        ranks_per_node = int(config.get('Run','ranks-per-node').partition('#')[0].strip())
    except:
        ranks_per_node = 32

    # ---------------- Start with problem-specific code ---------------------
    if problem is 'advection':
        # Supply a baseline dt and associated grid resolution so we can scale the
        # time step appropriately for different resolutions.
        dt_eff_res    = int(config.get('Problem','dt_eff_res').partition('#')[0].strip())
        dt_fixed      = float(config.get('Problem','dt_fixed').partition('#')[0].strip())

        # Replicated problem
        try:
            d = config.get('Run','duplicate').partition('#')[0].strip()
            duplicate  = d in ['T','True','1']
        except:
            duplicate = False

        # Uniform or adaptive
        use_fixed_dt = True
        outstyle = 3

        # Beta : Adapt proc count
        adapt_proc_count = False

    elif problem is 'shockbubble':
        initial_dt   = float(config.get('Run','initial_dt').partition('#')[0].strip())
        # smooth_level0  = maxlevel0-1   # Smooth at all levels
        smooth_level0    = float(config.get('Run','smooth-level').partition('#')[0].strip())
        try:
            nout0 = int(config.get('Run','nout').partition('#')[0].strip())
            tfinal0 = initial_dt*nout0
        except:
            nout0 = 1

        # Beta : adapt proc count
        try:
            num_grids_per_proc = int(config.get('Run','num_grids_per_proc').partition('#')[0].strip())
            adapt_proc_count = True
        except:
            adapt_proc_count = False

        use_fixed_dt = True
        outstyle = 3
        nstep0 = nout0
        duplicate = False

        # number of time steps
        try:
            coeff_C      = float(config.get('Problem','adapt_coeff_C').partition('#')[0].strip())
            coeff_D      = float(config.get('Problem','adapt_coeff_D').partition('#')[0].strip())
            # dt0  = float(config.get('Problem','dt_coarse_est').partition('#')[0].strip())
        except:
            pass

    elif problem is 'slotted_disk':
        initial_dt   = float(config.get('Run','initial_dt').partition('#')[0].strip())
        # smooth_level0  = maxlevel0-1   # Smooth at all levels
        smooth_level0    = float(config.get('Run','smooth-level').partition('#')[0].strip())
        nout0 = int(config.get('Run','nout').partition('#')[0].strip())
        nstep0 = int(config.get('Run','nstep').partition('#')[0].strip())
        tfinal0 = initial_dt*nout0

        # Beta : adapt proc count
        try:
            num_grids_per_proc = int(config.get('Run','num_grids_per_proc').partition('#')[0].strip())
            adapt_proc_count = True
        except:
            adapt_proc_count = False

        use_fixed_dt = True
        outstyle = 3
        duplicate = False

        # number of time steps
        try:
            coeff_A      = float(config.get('Problem','adapt_coeff_A').partition('#')[0].strip())
            coeff_B      = float(config.get('Problem','adapt_coeff_B').partition('#')[0].strip())
        except:
            pass

    # ---------------- Done with problem-specific code ---------------------

    try:
        verbosity = config.get('Run','verbosity').partition('#')[0].strip()
    except:
        verbosity = 'essential'

    # Mode : uniform or adaptive ?
    mode    = config.get('Run','mode').partition('#')[0].strip()


    # Subcycling (used only in adaptive case)
    sc = config.get('Run','subcycle').partition('#')[0].strip()
    subcycle = sc in ['T','True','1']

    try:
        nw = config.get('Run','weightedp').partition('#')[0].strip()
        weightedp = nw in ['T','True','1']
    except:
        weightedp = False

    try:
        gparea = config.get('Run','ghost_patch_pack_area').partition('#')[0].strip()
        gparea = gparea in ['T','True','1']
    except:
        gparea = True

    # Advance one step?
    try:
        ol = config.get('Run','advance-one-step').partition('#')[0].strip()
        one_level = ol in ['T','True','1']
    except:
        one_level = False

    # Use Maxlevel?
    try:
        ml = config.get('Run','outstyle-uses-maxlevel').partition('#')[0].strip()
        use_maxlevel = ml in ['T','True','1']
    except:
        use_maxlevel = False

    # these should eventually go away ...
    scaling = 'strong'   # Get rid of these eventually?
    scale_uniform = False

    # ----------------------------------------
    # Other inputs needed by the user
    # ----------------------------------------

    # ---------------- Start of problem-specific code ---------------------
    if problem is 'advection':
        # Figure out dt needed for first run in this series
        eff_res0 = mx0*2**minlevel0

        dt0 = dt_fixed/(float(eff_res0)/float(dt_eff_res))

        steps_coarse0 = tfinal0/dt0

        if use_maxlevel:
            steps_coarse = steps_coarse0*2**(maxlevel0-minlevel0)
        else:
            steps_coarse = steps_coarse0

        tol = 3e-15
        if abs(steps_coarse - np.round(steps_coarse0)) > tol:
            print("steps_coarse is not an integer; nout = %12.4e\n" % \
                (steps_coarse0-np.round(steps_coarse0)))
            sys.exit()
            steps_coarse = int(steps_coarse0)

        nout0 = steps_coarse0  # Used for outstyle == 3
        nstep0 = nout0  # Only output final time step.
        initial_dt = dt0

    elif problem is 'shockbubble':
        smooth_level0  = 0
        steps_coarse0 = nout0
        dt0 = tfinal0/steps_coarse0   # coarse grid time step
    elif problem is 'slotted_disk':
        steps_coarse0 = nout0
        dt0 = initial_dt
        pass

    # ---------------- Done with problem-specific code ---------------------


    # ----------------------------------------
    # Everything should be automatic from here
    # Four types of runs ;
    #   -- uniform with weak scaling
    #   -- uniform with strong scaling
    #   -- adaptive with weak scaling
    #   -- adaptive with strong scaling
    # ----------------------------------------

    R = np.arange(0,njobs)

    if mode == 'adapt':
        adapt_factor = np.exp(coeff_B*maxlevel0 + coeff_A)
    else:
        adapt_factor = np.ones(R.shape)

    mx       = np.empty(R.shape)
    nout     = np.empty(R.shape)
    nstep    = np.empty(R.shape)
    minlevel = np.empty(R.shape)
    maxlevel = np.empty(R.shape)
    tfinal   = np.empty(R.shape)
    dt       = np.empty(R.shape)
    mi       = np.empty(R.shape)
    mj       = np.empty(R.shape)
    smooth_level  = np.empty(R.shape)
    steps_coarse  = np.empty(R.shape)

    # These are fixed for all possible runs.
    mx.fill(mx0)
    nout.fill(nout0)   # set to 1 in shockbubble problem
    nstep.fill(nstep0)
    steps_coarse.fill(steps_coarse0)

    # minlevel is fixed
    minlevel.fill(minlevel0)
    maxlevel.fill(maxlevel0)

    mi.fill(mi0)
    mj.fill(mj0)

    # Run entire simulation, to get true measure of adaptivity
    tfinal.fill(tfinal0)
    dt.fill(dt0)    # Coarse grid dt; is scaled to finer levels.

    if problem is 'advection':
        smooth_level = maxlevel-1

    if problem is 'shockbubble':
        smooth_level.fill(smooth_level0)

    if problem is 'slotted_disk':
        smooth_level = maxlevel-1


    nout_uniform = steps_coarse*2**(maxlevel-minlevel)  # Number of fine grid steps
    num_grids_total = mi0*mj0*(4**maxlevel)  # total on all procs, all blocks

    if adapt_proc_count:
        # Doesn't work yet
        proc0 = np.ceil(np.exp(coeff_B*maxlevel0 + coeff_A)*num_grids_total/num_grids_per_proc)
        f = 4*np.exp(coeff_B)
        procs = proc0*f**R
    else:
        procs = proc0*4**R     # Key idea : Run on different processors

    if not duplicate:
        # grids_per_proc = adapt_factor*mi0*mj0*(4**maxlevel)/procs
        grids_per_proc = adapt_factor*num_grids_total/procs
    else:
        grids_per_proc = adapt_factor*(4**maxlevel)/(procs/procs[0])

    num_advances = adapt_factor*nout_uniform*num_grids_total
    t = (num_advances*time_per_grid)/procs  # Assuming ideal scaling

    # eff_res = mx*(2**maxlevel)*mi0
    eff_res = np.round(np.sqrt((mx*2**maxlevel)**2*mi0*mj0))

    # ------------------------------------------
    # Start creating files.
    # ------------------------------------------
    fmt_str_numeric = "# %6d %4d %7d %7d %6d %6d %12.8g %12.8g %8d %12d %8.1f"
    fmt_str_header  = "# %6s %4s %7s %7s %6s %6s %12s %12s %8s %12s %8s"
    tuple_str = ('p','mx','minlev','maxlev','nout','nstep', 'dt','tfinal',
                 'eff_res','grids/proc','t')
    header_line = "# " + "-"*99

    # Output to console
    print(header_line)
    print(fmt_str_header % tuple_str)
    print(header_line)

    # Output to job files (jobs.sh)
    jobfile = open('jobs.sh','w')
    jobfile.write("#!/usr/bin/env bash\n")
    jobfile.write("#\n")
    jobfile.write("# Time per grid : %8.2e\n" % time_per_grid)
    jobfile.write("#\n")
    jobfile.write(header_line)
    jobfile.write("\n")
    jobfile.write((fmt_str_header + "\n") % tuple_str)
    jobfile.write(header_line)
    jobfile.write("\n")

    filein = open('results.in','w')
    filein.write(header_line)
    filein.write("\n")
    filein.write((fmt_str_header + "\n") % tuple_str)
    filein.write(header_line)
    filein.write("\n")

    cpu_hours = 0
    for i,p in enumerate(procs):
        if grids_per_proc[i] < 4:
            print("Too few grids per proc. No more output files will be created.")
            sys.exit()

        if duplicate:
            level_inc = (np.log(float(mi0))/np.log(float(2.0))).astype(int)
        else:
            level_inc = 0

        prt_tuple = (p,mx[i],minlevel[i]+level_inc,maxlevel[i]+level_inc,
                     steps_coarse[i],nstep[i],dt[i],
                     tfinal[i],eff_res[i], grids_per_proc[i],t[i])

        # Print to console and jobs file
        print(fmt_str_numeric % prt_tuple)
        jobfile.write((fmt_str_numeric + "\n") % prt_tuple)
        # Skip comments
        filein.write(("  " + fmt_str_numeric[2:] + "\n") % prt_tuple)

        if scheduler == 'll':
            cpu_hours = cpu_hours + np.max([p,32])*t[i]/(3600.0)
        else:
            cpu_hours = cpu_hours + t[i]/3600.0

        inifile = "ex_%05d.ini" % p
        ini_file = open(inifile,'w')
        ini_file.write("[user]\n")
        ini_file.write("     example = %d      # mapping\n" % example)
        ini_file.write("\n")
        ini_file.write("[Options]\n")
        ini_file.write("# Grid dimensions\n")
        ini_file.write("    mx = %d\n" % mx[i])
        ini_file.write("    my = %d\n" % mx[i])
        ini_file.write("\n")
        ini_file.write("    minlevel = %d\n" % minlevel[i])
        ini_file.write("    maxlevel = %d\n" % maxlevel[i])

        ini_file.write("    regrid_interval = %d\n" % regrid_interval)
        ini_file.write("    smooth-level = %d\n" % (smooth_level[i]))

        ini_file.write("\n")
        ini_file.write("    tfinal = %16.8e\n" % tfinal[i])
        if use_fixed_dt:
            ini_file.write("    use_fixed_dt = T\n")
        else:
            ini_file.write("    use_fixed_dt = F\n")

        ini_file.write("    initial_dt = %20.16e\n" % initial_dt)
        ini_file.write("\n")
        ini_file.write("    outstyle = %d\n" % outstyle)
        ini_file.write("    nout = %d\n" % nout[i])
        ini_file.write("    nstep = %d\n" % nstep[i])  # Not used if outstyle == 1
        if mode == 'adapt':
            if subcycle:
                ini_file.write("    subcycle = T\n")
            else:
                ini_file.write("    subcycle = F\n")

            if one_level:
                ini_file.write("    advance-one-level = T\n");
            else:
                ini_file.write("    advance-one-step = F\n")

            if use_maxlevel:
                ini_file.write("    outstyle-uses-maxlevel = T\n");
            else:
                ini_file.write("    outstyle-uses-maxlevel = F\n")

        else:
            ini_file.write("    subcycle = F\n")


        ini_file.write("    init_ghostcell = F\n")
        ini_file.write("\n")

        ini_file.write("# Subcycling\n");
        if weightedp:
            ini_file.write("    weighted_partition = T\n")
        else:
            ini_file.write("    weighted_partition = F\n")


        if gparea:
            ini_file.write("    ghost_patch_pack_area = T\n")
        else:
            ini_file.write("    ghost_patch_pack_area = F\n")

        ini_file.write("\n")

        # Other things which should not be set of timing runs
        ini_file.write("# File and console IO\n")
        ini_file.write("    verbosity = %s\n" % verbosity)
        ini_file.write("    serialout = F\n")
        ini_file.write("    vtkout = 0\n")
        ini_file.write("    tikzout = F\n")
        ini_file.write("\n")

        ini_file.write("# Debugging and diagnostics\n")
        ini_file.write("    conservation-check = F\n")
        ini_file.write("    run-user-diagnostics = F\n")
        ini_file.write("    trapfpe = F\n")
        ini_file.write("    mpi_debug = F\n")
        ini_file.write("\n")

        ini_file.write("# Mapping options\n")
        ini_file.write("    mi = %d\n" % (mi[i]))
        ini_file.write("    mj = %d\n" % (mj[i]))


        ini_file.close()

        if scheduler == 'll':
            # JUQUEEN (BlueGene/Q)
            procfile = "p_%05d.ll" % p
            proc_file = open(procfile,'w')
            proc_file.write("#@ job_name = ex_%05d.ini\n" % p)
            proc_file.write("\n")
            proc_file.write("#@ comment = \"%s example : eff. resolution = %d x %d\"\n" % \
                            (execname.capitalize(),eff_res[i],eff_res[i]))
            proc_file.write("#@ error = %s_%05d.o$(jobid)\n" % (execname,p))
            if problem == 'slotted_disk':
                proc_file.write("#@ output = %s_%05d.o$(jobid)\n" % ('slotted-disk',p))
            else:
                proc_file.write("#@ output = %s_%05d.o$(jobid)\n" % (execname,p))

            proc_file.write("\n")
            proc_file.write("#@ environment = COPY_ALL\n")
            trun = np.min([60**2,1.1*t[i]])  # Keep jobs under 1 hour
            if p <= 1024:
                if trun > 60*30:
                    # print "Time on proc %d may exceed allowable time of 30 minutes" % p
                    pass
                proc_file.write("#@ wall_clock_limit = 00:30:00\n")
            else:
                if trun > 60**2:
                    # print "Time may exceed 1 hour"
                    pass
                proc_file.write("#@ wall_clock_limit = 01:00:00\n")

            proc_file.write("#@ notification = error\n")
            proc_file.write("#@ notify_user = donnacalhoun@boisestate.edu\n")
            proc_file.write("#@ job_type = bluegene\n")
            rpn = np.min([ranks_per_node,p])
            bgsize = np.max([32,np.ceil(p/rpn)])
            proc_file.write("#@ bg_size = %d\n" % bgsize)
            proc_file.write("#@ queue\n")
            proc_file.write("\n")
            proc_file.write(("runjob --ranks-per-node %d --np %d : " \
                             + "/homec/hbn26/hbn263/projects/forestclaw-build/local/bin/%s " \
                             + "--inifile=ex_%05d.ini\n") %
                            (np.min([rpn,p]),p,execname,p))
            proc_file.close()

        elif scheduler == 'pbs':
            # Kestrel (Linux Cluster)
            procfile = "p_%05d.pbs" % p
            proc_file = open(procfile,'w')
            proc_file.write("#!/bin/bash\n")
            trun = np.min([24.0*60**2,2*t[i]])  # Keep jobs under 24 hours
            if trun < 60*30:
                proc_file.write("#PBS -l walltime=00:30:00\n")
            else:
                h = trun//3600
                m = (trun-3600*h)//60
                s = trun - 3600*h - 60*m
                proc_file.write("#PBS -l walltime=%02d:%02d:%02d\n" % (h,m,s))

            if p <= 16:
                proc_file.write("#PBS -l select=1:ncpus=%d:mpiprocs=%d:mem=2gb\n" % (p,p))
            elif p == 64:
                proc_file.write("#PBS -l select=4:ncpus=16:mpiprocs=16:mem=2gb\n")
            else:
                print("PBS scheduler : Unexpected number of processors; p = %d\n" % p)
                sys.exit()

            proc_file.write("\n")
            proc_file.write("#PBS -N %s_%05d\n" % (execname,p))
            proc_file.write("#PBS -q batch\n")
            proc_file.write("#PBS -j oe\n")
            proc_file.write("#PBS -r n\n")
            proc_file.write("#PBS -l place=scatter:excl\n")
            proc_file.write("\n")
            proc_file.write("module load shared\n")
            proc_file.write("module load pbspro\n")
            proc_file.write("module load gcc/4.8.1\n")
            proc_file.write("module load openmpi/1.8.5\n")
            proc_file.write("\n")
            proc_file.write("cd $PBS_O_WORKDIR\n")
            exec_dir = "/home/donnacalhoun/projects/forestclaw-build/local/bin"
            proc_file.write("mpirun %s/torus --inifile=%s\n" % (exec_dir,inifile))

        elif scheduler == 'osx':
            # Mac (laptop or desktop). Write a Python script that calls
            # mpirun as a subprocess.
            procfile = "p_%05d.py" % p
            proc_file = open(procfile,'w')
            proc_file.write("# comment = \"%s example : eff. resolution = %d x %d\"\n" % \
                            (execname.capitalize(),eff_res[i],eff_res[i]))

            proc_file.write("import sys\n")
            proc_file.write("import os\n");
            proc_file.write("import subprocess\n")
            proc_file.write("import random\n")
            proc_file.write("\n")

            arg_list = ["mpirun", "-n", "$p", "$execname","--inifile=$inifile"]
            sep = "\",\""
            s = Template("arg_list = [\"" + sep.join(arg_list) + "\"]\n")
            alist = s.substitute(p=str(p),execname=execname,inifile=inifile)

            proc_file.write(alist)
            proc_file.write("jobid = random.randint(1000,9999)\n")
            ostr = Template("outfile = \"$execname.o%d\" % (jobid)\n")
            ostr = ostr.substitute(execname="%s_%05d" % (execname,p))
            proc_file.write(ostr)
            proc_file.write("f = open(outfile,'w')\n")
            proc_file.write("po = subprocess.Popen(arg_list,stdout=f)\n")
            pstr = "print(\"Starting process %d with jobid %d on $p " + \
                    "processor(s).\" % (po.pid,jobid))\n"
            s = Template(pstr)
            pstr = s.substitute(p=str(p))
            proc_file.write(pstr)
            proc_file.write("po.wait()\n")
            proc_file.close()


    print(header_line)
    print("Estimated core-h (hours) : %6.2f (%4.2f%%)" % (cpu_hours,100*cpu_hours/350000.0))


    jobfile.write(header_line)
    jobfile.write("\n")
    jobfile.write("# Estimated core-h (hours) : %6.2f (%4.2f%%)" % \
                  (cpu_hours,100*cpu_hours/3500000.0))
    jobfile.write("\n\n");

    filein.write(header_line)
    filein.write("\n")
    filein.write("# Estimated core-h (hours) : %6.2f (%4.2f%%)" % \
                  (cpu_hours,100*cpu_hours/3500000.0))
    filein.write("\n\n")
    filein.close()

    for i,p in enumerate(procs):

        if scheduler == 'll':
            fname = "p_%05d.ll" % p
            jobfile.write("llsubmit %s\n" % fname)
        elif scheduler == 'pbs':
            fname = "p_%05d.pbs" % p
            jobfile.write("qsub %s\n" % fname)
            pass
        elif scheduler == 'osx':
            fname = "p_%05d.py" % p
            jobfile.write("python %s\n" % fname)

    jobfile.close()

    print("")
    print("Configuration files written : ")
    for i,p in enumerate(procs):
        inifile = "ex_%05d.ini" % p
        print("    %s" % inifile)

    print("")
    print("jobs submitted from jobs.sh")
    for i,p in enumerate(procs):
        if scheduler == 'll':
            procfile = "p_%05d.ll" % p
        elif scheduler == 'pbs':
            procfile = "p_%05d.pbs" % p
        elif scheduler == 'osx':
            procfile = "p_%05d.py" % p
        print("    %s" % procfile)


def launch_jobs(N=1):

    import subprocess
    for i in range(N):
        po = subprocess.call(['bash','jobs.sh'])

