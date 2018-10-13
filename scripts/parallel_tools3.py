import os
import sys
import numpy as np
from string import Template
import re
import configparser
from pdb import set_trace


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

        # Replicated problem
        try:
            replicate_factor = int(config.get('Run','replicate-factor').partition('#')[0].strip())
        except:
            replicate_factor = 1

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

    R = np.arange(njobs)

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
    if duplicate:
        rf = replicate_factor
        num_grids_total = rf*rf*(4**maxlevel)  # total on all procs, all blocks

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

        # Levels are incremented above
        if duplicate:
            level_inc = 0
            rf = replicate_factor
            while rf != 1:
                rf /= 2;
                level_inc +=1
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
        if duplicate:
            ini_file.write("     minlevel-base = {:d}\n".format(minlevel0))
            ini_file.write("     maxlevel-base = {:d}\n".format(maxlevel0))
            ini_file.write("     replicate-factor = {:d}\n".format(replicate_factor))

        ini_file.write("\n")
        ini_file.write("[clawpatch]\n")
        ini_file.write("# Grid dimensions\n")
        ini_file.write("    mx = %d\n" % mx[i])
        ini_file.write("    my = %d\n" % mx[i])
        ini_file.write("\n")
        ini_file.write("[Options]\n")
        if not duplicate:
            ini_file.write("    minlevel = {:d}\n".format(minlevel[i]))
            ini_file.write("    maxlevel = {:d}\n".format(maxlevel[i]))

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
        ini_file.write("    output = F\n")
        ini_file.write("\n")

        ini_file.write("# Debugging and diagnostics\n")
        ini_file.write("    conservation-check = F\n")
        ini_file.write("    run-user-diagnostics = F\n")
        ini_file.write("    trapfpe = F\n")
        ini_file.write("    mpi_debug = F\n")
        ini_file.write("\n")

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
                         'ADVANCE_B4STEP2',
                         'ADVANCE_STEP2',
                         'CUDA_MEMCOPY_H2H',
                         'CUDA_MEMCOPY_H2D',
                         'CUDA_MEMCOPY_D2H',
                         'OUTPUT']

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
                  'interp','average','adv_b4step2','adv_step2','memcopy_h2h',
                  'memcopy_h2d','memcopy_d2h','output')
    float_width = 12*len(float_list)
    float_str = ("{:>12s}"*len(float_list)).format(*float_list)


    # Compile data to this file
    line_len = option_width + int_width + float_width
    header_str = option_str_header + int_str + float_str

    # Write out header to file
    resultsfile = open(results_file,'w')
    # resultsfile.write("# " + "-"*line_len)
    # resultsfile.write("\n")

    resultsfile.write("# " + header_str + "\n")

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


# Read in an array of results that can be used for both weak
# and strong scaling.
def read_results_files(dir_list, subdir = None, results_in = None,
                       results_file='results.out',execname=None):

    import re

    dirs = []
    for d in dir_list:
        dirs.append("run_%03d" % d)

    read_input = not (results_in == None)

    if read_input:
        # Get data from results.in only
        pattern = re.compile(results_in)
        results_file = results_in
    else:
        if not execname == None:
            pattern = re.compile("%s_[0-9]{5}\.o[0-9]*" % execname)
        else:
            pattern = re.compile(".*_[0-9]{5}\.o[0-9]*")

    mx1 = []
    procs1 = []
    levels1 = []
    t = []
    for results_dir in dirs:
        # Match file names like run_004
        if subdir == None:
            d = results_dir
        else:
            d = os.path.join(results_dir,subdir)


        files = os.listdir(d)
        for f in files:
            if re.match(pattern,f):
                # Load data file to read in (mx,levels) for this file
                if subdir == None:
                    rf = os.path.join(results_dir,results_file)
                else:
                    rf = os.path.join(results_dir,subdir,results_file)

                if not os.path.exists(rf):
                    print("File %s not found in %s" % (results_file,f))
                    continue

                data = np.loadtxt(rf)
                if data.ndim == 1:
                    # Needed in case we only have 1 row in results file.
                    data = np.reshape(data,(-1,len(data)))

                if read_input:
                    procs1.extend(data[:,0])
                    mx1.extend(data[:,1])
                    levels1.extend(data[:,3])
                else:
                    procs1.extend(data[:,1])    # Must be updated if results file is changed
                    mx1.extend(data[:,2])       # Check below
                    levels1.extend(data[:,6])   # Check  below

    t = zip(mx1,procs1,levels1)  # All possible combos

    # Create array of all (procs/levels).  Not every combo will have data.
    mx = sorted(list(set(mx1)))
    procs = sorted(list(set([x[1] for x in t])))
    levels = sorted(list(set([x[2] for x in t])))

    jobs = dict.fromkeys(mx)
    for m in mx:
        jobs[m] = dict.fromkeys(procs)
        for p in procs:
            jobs[m][p] = dict.fromkeys(levels)


    # Initialize dictionaries for all values
    for results_dir in dirs:
        rundir = int(results_dir.partition('_')[2])
        if subdir == None:
            d = results_dir
        else:
            d = os.path.join(results_dir,subdir)

        files = os.listdir(d)
        for f in files:
            if re.match(pattern,f):

                # Load data file to read in (mx,levels) for this file
                if subdir == None:
                    rf = os.path.join(results_dir,results_file)
                else:
                    rf = os.path.join(results_dir,subdir,results_file)

                if not os.path.exists(rf):
                    print("File %s not found in %s" % (results_file,f))
                    continue

                data = np.loadtxt(rf)

                if data.ndim == 1:
                    data = np.reshape(data,(-1,len(data)))

                for row in data:
                    if read_input:
                        p = row[0]
                        m = row[1]
                        l = row[3]

                        s = int(results_dir.partition('_')[2])

                        job = {}
                        job["rundir"]         = rundir
                        job["procs"]          = row[0]
                        job["mx"]             = row[1]
                        job["minlevel"]       = row[2]
                        job["maxlevel"]       = row[3]
                        job["nout"]           = row[4]
                        job["nstep"]          = row[5]
                        job["dt"]             = row[6]
                        job["tfinal"]         = row[7]
                        job["effres"]         = row[8]
                        job["grids_per_proc"] = row[9]
                        job["walltime"]       = row[10]

                    else:
                        p = row[1]
                        m = row[2]    # Index of mx in data array
                        l = row[6]    # Index of maxlevel

                        if len(row) < 32:
                            print("read_results_file (parallel_tools.py)")
                            print("Length of row from %s is incorrect" %(rf))
                            print("It should be 32 instead of %d" % (len(row)))
                            sys.exit()


                        job = {}
                        job["rundir"]          = rundir
                        job["jobid"]           = row[0]
                        job["procs"]           = row[1]
                        job["mx"]              = row[2]   # Used above in computing m
                        job["mi"]              = row[3]
                        job["mj"]              = row[4]
                        job["minlevel"]        = row[5]
                        job["maxlevel"]        = row[6]   # Used above in computing l
                        job["nout"]            = row[7]
                        job["tfinal"]          = row[8]
                        job["advance_steps"]   = row[9]
                        job["grids_per_proc"]  = row[10]
                        if job["grids_per_proc"] == 0:    # avoid log(0)
                            job["grids_per_proc"] = 0.1

                        job["local_boundary"]  = row[11]
                        job["remote_boundary"] = row[12]
                        job["interior"]        = row[13]
                        job["init"]            = row[14]
                        job["advance"]         = row[15]
                        job["ghostfill"]       = row[16]
                        job["regrid"]          = row[17]
                        job["ghostpatch_comm"] = row[18]
                        job["adapt_comm"]      = row[19]
                        job["partition_comm"]  = row[20]
                        job["cfl_comm"]        = row[21]
                        job["walltime"]        = row[22]
                        job["ghostfill_copy"]  = row[23]
                        job["local"]           = row[24]
                        job["comm"]            = row[25]
                        job["partition"]       = row[26]
                        job["ghostpatch_build"]= row[27]
                        job["step1"]           = row[28]
                        job["step2"]           = row[29]
                        job["step3"]           = row[30]
                        job["copy"]            = row[31]
                        job["interp"]          = row[32]
                        job["average"]         = row[33]
                        job["local_ratio"]     = row[34]
                        job["remote_ratio"]    = row[35]

                    jobs[m][p][l] = job

    return jobs


def print_jobs(jobs,val2plot,fmt_int=False):

    mx = sorted(jobs.keys())

    int_list = ['jobid','procs','mx','minlevel','maxlevel','nout', 'advance_steps',
                'rundir','mi','mj']

    for m in mx:
        procs = sorted(jobs[m].keys())
        levels = []
        for p in procs:
            levels += jobs[m][p]

        levels = np.array(sorted(set(levels))).astype(int)
        data = np.empty((len(procs),len(levels)))

        rundir_data = np.empty(levels.shape)

        for i,p in enumerate(procs):
            for j,l in enumerate(levels):
                try:
                    job = jobs[m][p][l]
                except:
                    data[i,j] = None
                    continue                

                if job == None:
                    # This (m,p,l) combo doesn't exist;  p
                    data[i,j] = np.nan
                    continue
                else:
                    rundir_data[j] = job["rundir"]

                if isinstance(val2plot,str):
                    try:
                        d = job[val2plot]
                    except:
                        print("job[\"%s\"] not a valid variable" % (val2plot))
                        pdb.set_trace()
                        sys.exit()
                else:
                    import types
                    if isinstance(val2plot,types.FunctionType):
                        try:
                            d,fmt_int = val2plot(job,mx=m,proc=p,level=l,all=jobs)
                        except:
                            print("print_jobs : Problem with function %s" % (val2plot.__name__))
                            sys.exit()


                try:
                    data[i,j] = np.average(d)
                except:
                    print("job (mx=%d,proc=%d,level=%d) is empty" % (m,p,l))
                    data[i,j] = np.nan

        # Create an array that interleaves levels and rundirs so that
        # we can print out both in the header, e.g.
        # 2 (10)   3 (11)    4 (12)
        header_data = zip(levels,rundir_data.astype(int))

        # Header
        print("")
        print(val2plot)
        linelen = 8 + 3 + 13*len(levels)
        line = "-"*linelen
        print(line)
        mxstr = "mx = %d" % m

        # header_format = "{mxstr:>8}" + "{sep:>3}" + "{:>12d}" * len(levels)
        # print header_format.format(mxstr = mxstr, sep="|",*header_data)

        s = '{mxstr:>8}'.format(mxstr=mxstr) + '{sep:>3}'.format(sep='|')
        for t in header_data:
            s = s +  '{level:>7d} ({rundir:3d})'.format(level=t[0],rundir=t[1])
        print(s)
        print(line)
        # Format integer and float values separately
        if (val2plot in int_list) or fmt_int:
            row_format = "{:>8d}" + "{sep:>3}" + "{:>13.0f}"*len(levels)
        else:
            row_format = "{:>8d}" + "{sep:>3}" + "{:>13.2e}"*len(levels)
        for p, d in zip(procs, data):
            print(row_format.format(int(p), sep="|", *d))

        print(line)


def plot_results(jobs,start_point,val2plot='walltime',
                 scaling='strong',scale_uniform=False,
                 ideal_slope = True,efficiency=False,
                 scatter=False,model='normal'):

    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.lines import Line2D

    import seaborn as sns

    ax = plt.gca()
    # Or we can set the colormap explicitly.
    cmap = sns.color_palette("Set1",10)
    cm = [cmap[i] for i in [0,5,1,4,2,7,3]]
    from cycler import cycler
    ax.set_prop_cycle(cycler('color',cm))


    markers = {}
    markers["walltime"] = r'$\clubsuit$'
    markers["advance"] = u's'
    markers["ghostfill"] = u'h'
    markers["regrid"] = u'v'
    markers["ghostpatch_comm"] = u'^'
    markers["advance_steps"] = u'o'
    markers["cfl_comm"] = u'p'

    # colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')*2

    mx = sorted(jobs.keys())
    ph = []
    mb_data = []

    # each start point generates a plot, either down a column (strong), across a
    # diagonal (weak) or from back to front (block size comparison)
    tlist = []
    gpplist = []
    ylist = []
    proclist = []
    for i,sp in enumerate(start_point):
        m = sp[0]
        p = sp[1]
        l = sp[2]
        try:
            job = jobs[m][p][l]
        except:
            print("Job entry not found;  m = %d; p = %d; l = %d" % (m,p,l))            
            import pdb
            pdb.set_trace()
            job = None


        # "headers" for this mx table
        procs = sorted(jobs[m].keys())
        # levels = sorted(jobs[m][p].keys())
        # levels = np.array(sorted(set(levels))).astype(int)
        levels = []
        for p in procs:
            levels += jobs[m][p].keys()

        levels = np.array(sorted(set(levels))).astype(int)

        if scaling  == 'weak':  # Work down/right across diagonal
            mx_sub = [sp[0] for p in procs if p >= sp[1]]
            procs_sub =  [p for p in procs if p >= sp[1]]
            # procs_sub = procs
            levels_sub = [l for l in levels if l >= sp[2]]

        elif scaling == 'strong':  # Work down colum
            # assume that mx is fixed
            mx_sub = [sp[0] for p in procs if p >= sp[1]]
            procs_sub = []
            for p in procs:
                try:
                    job = jobs[m][p][l]
                except:
                    pass
                else:
                    procs_sub += [p]
            # procs_sub = [p for p in procs if p >= sp[1]]
            levels_sub = [sp[2] for p in procs_sub]


        elif scaling == 'resolution':
            # scaling as resolution is increased
            levels_sub = [l for l in levels if l >= sp[2]]
            mx_sub = [sp[0] for l in levels_sub]
            procs_sub = [sp[1] for l in levels_sub]

        t = list(zip(mx_sub,procs_sub,levels_sub))

        # Internal plotting routine
        phandle,mb,gpp1,y_avg1,procs1 = plot_results_internal(val2plot,t,
                                                              jobs,scaling,markers,
                                                              scale_uniform,efficiency,
                                                              ideal_slope,
                                                              scatter)

        # phandle.set_color(colors[i])
        plt.draw()
        ph.append(phandle)
        mb_data.append(mb)

        tlist.extend(t)
        gpplist.extend(gpp1)
        ylist.extend(y_avg1)
        proclist.extend(procs1)

    # Set up axis
    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%d'))

    A = None
    B = None
    alpha = None

    if scatter:
        plt.xscale('log')
        plt.xlim([2**(-0.5),1e4])
        xticks = ticker.MultipleLocator(base=10)
        xticks = ticker.MaxNLocator(nbins=7)
        ax.xaxis.set_minor_locator(xticks)

        xformatter = ticker.FormatStrFormatter('%d')
        ax.xaxis.set_major_formatter(xformatter)
        ax.set_xlabel("Grids per processor",fontsize=16)

        yformatter = ticker.FormatStrFormatter('%5.2e')
        ax.yaxis.set_major_formatter(yformatter)
        ax.set_ylabel("Cost per grid",fontsize=16)

        if len(gpplist) > 1:
            ranks_per_node = -1    # Set to -1 to get single curve through all data
            gpp = np.array(gpplist)
            y = np.array(ylist)
            z = zip(gpp,y,proclist)
            gpp      = np.array([x[0] for x in z if x[0] > 1]).astype(int)
            y        = np.array([x[1] for x in z if x[0] > 1])
            proclist = np.array([x[2] for x in z if x[0] > 1]).astype(int)

            z = zip(gpp,y,proclist)
            gpp1 = np.array([t[0] for t in z if t[2] <= ranks_per_node])
            y1   = np.array([t[1] for t in z if t[2] <= ranks_per_node])

            if len(gpp1) > 1:
                xs = np.logspace(np.log10(np.min(gpp1)),np.log10(np.max(gpp1)),50)
                ys,A,alpha,B = scatter_fit(gpp1,y1,xs,model)
                plt.plot(xs,ys,'k--')

            gpp2 = np.array([t[0] for t in z if t[2] > ranks_per_node])
            y2   = np.array([t[1] for t in z if t[2] > ranks_per_node])
            if len(gpp2) > 1:
                xs = np.logspace(np.log10(np.min(gpp2)),np.log10(np.max(gpp2)),50)
                ys,A,alpha,B = scatter_fit(gpp2,y2,xs,model)
                plt.plot(xs,ys,'k-')

            plt.draw()

    else:
        if scaling == 'resolution':
            levels = list(set([x[2] for x in tlist]))
            ax.xaxis.set_major_locator(plt.FixedLocator(levels))
            ax.set_xlabel("Levels",fontsize=16)
            l1 = np.min(levels)
            l2 = np.max(levels)
            plt.xlim([l1-0.5,l2+0.5])

        else:
            ax.xaxis.set_major_locator(plt.FixedLocator(procs))
            ax.set_xlabel("Processor count",fontsize=14)
            p1 = np.log(np.min(procs))/np.log(4)
            p2 = np.log(np.max(procs))/np.log(4)
            plt.xlim([4**(p1-0.5), 4**(p2+0.5)])

        if efficiency:
            ax.set_yscale('linear')
            ax.set_yticks(range(0,110,10),minor=True)
            plt.grid(b=True,which='major')
            plt.grid(b=True,which='minor',axis='y')
            # plt.minorticks_on()
            plt.ylim([0,110])

    plt.grid(b=True,which='minor',axis='y')
    plt.grid(b=True,which='major',axis='x')
    plt.setp(ax.get_xticklabels(),fontsize=11)
    plt.setp(ax.get_yticklabels(),fontsize=11)

    plt.draw()

    return ph,mb,t,A,alpha,B


def plot_results_internal(val2plot,t,jobs,scaling,markers,
                          scale_uniform=False, efficiency=False,
                          ideal_slope=None,scatter=False):

    import matplotlib.pyplot as plt

    if isinstance(val2plot,list):
        print("Only one val2plot allowed per plot")
        sys.exit()

    v = val2plot

    # mode = 'uniform'
    ph = []
    mx = sorted(jobs.keys())
    y_avg = []
    gpp = []

    # First collect y values (which may need to be averaged)
    for idx_tuple in t:

        m = idx_tuple[0]
        p = int(idx_tuple[1])
        l = idx_tuple[2]
        job = jobs[m][p][l]
        if job == None:
            y_avg = np.append(y_avg,np.nan)
            gpp = np.append(gpp,np.nan)
            continue

        if isinstance(v,str):
            try:
                data = np.average(job[v])  # Average multiple runs
            except:
                print("job[\"%s\"] not a valid variable" % (val2plot))
                sys.exit()
        else:
            try:
                # Call function val2plot
                val,_ = v(job,mx=m,proc=p,level=l,all=jobs)
                data = np.average(val)
            except:
                print("plot_results_internal : Problem with function '%s'" % (val2plot.__name__))
                sys.exit()

        y_avg = np.append(y_avg,data)
        gpp = np.append(gpp,job["grids_per_proc"])

    if efficiency:
        y_avg = 100*y_avg/y_avg[0]

    # Plot ideal scaling
    tp = list(zip([x[0] for x in t],
             [x[1] for x in t],
             [x[2] for x in t],y_avg,gpp))
    mx      = [x[0] for x in tp if not np.isnan(x[3])]
    procs1  = [x[1] for x in tp if not np.isnan(x[3])]
    levels1 = [x[2] for x in tp if not np.isnan(x[3])]
    y_avg1  = [x[3] for x in tp if not np.isnan(x[3])]
    gpp1  = [x[4] for x in tp if not np.isnan(x[4])]

    mb = None

    if ideal_slope is not None:
        try:
            if ideal_slope == 0:
                s = 'k--'
            else:
                s = 'k.--'

            R = np.array([np.log(float(x))/np.log(4.0) for x in procs1])
            if R[0] > 0:
                R = R - R[0]
            ideal = y_avg1[0]*2**(ideal_slope*R)
            if scaling is not  'resolution':
                plt.loglog(procs1,ideal,s,markersize=15)
            else:
                plt.semilogy(levels1,ideal,s,markersize=15)
        except:
            print("plot_results_internal: Could not plot ideal curve")

    try:
        mv = markers[v]
    except:
        mv = u'o'

    if not scatter:
        # Here is where we finally plot the data
        if scaling in ['strong','weak']:
            ph = plt.loglog(procs1,y_avg1,marker=mv,markersize=10)
        elif scaling in ['resolution']:
            ph = plt.semilogy(levels1,y_avg1,marker=mv,markersize=10)

    else:
        # Plot scatter of data
        ph = plt.plot(gpp1,y_avg1,linestyle='',marker=mv,markersize=10,
                      color='b')


    return ph[0],mb,gpp1,y_avg1,procs1







def run_results(execname,run_dir,duplicate):

    dirs = []
    for d in run_dir:
        dirs.append("run_%03d" % d)

    cwd = os.getcwd()

    for results_dir in dirs:
        d = results_dir
        os.chdir(d)
        print("Compiling results in %s" % d)
        compile_results(duplicate=duplicate)
        os.chdir(cwd)






