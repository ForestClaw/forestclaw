import os
import sys
import numpy as np
from string import Template
import re
import ConfigParser

def compare_runs(execname,rerun,iframe):
    import random
    import subprocess
    import shutil
    # Compare serial and parallel runs

    files = os.listdir(os.getcwd())
    pattern = re.compile("p_[0-9]{5}.py")
    proc_list = []
    for f in files:
        if re.match(pattern,f):
            # Do a run for each processor file
            proc_list.append(int(f.partition('_')[2].partition('.py')[0]))

    if rerun:
        for p in proc_list:
            arg_list = ["mpirun","-n","%d" % p,execname,"--inifile=ex_%05d.ini" % p,
                        "--serialout=T"]

            jobid = random.randint(1000,9999)
            outfile = "%s_%05d.o%d" % (execname,p,jobid)
            f = open(outfile,'w')
            po = subprocess.Popen(arg_list,stdout=f)
            print "Starting process %d with jobid %d on %d processor(s)." % (po.pid,jobid,p)
            po.wait()

            rundir = 'run_%05d' % p
            shutil.rmtree(rundir,True)
            os.mkdir(rundir)
            pattern = re.compile("fort.[tq][0-9]{4}")

            fortfiles = os.listdir(os.getcwd())
            for f in fortfiles:
                if re.match(pattern,f):
                    shutil.move(f,rundir)

    else:
        for i,p in enumerate(proc_list):
            for q in proc_list[i+1:]:
                compare_dir = "compare_%05d_%05d" % (p,q)
                try:
                    os.mkdir(compare_dir)
                except:
                    if rerun:
                        shutil.rmtree(compare_dir,True)
                        print "Creating directory %s" % compare_dir
                else:
                    print "Using directory %s" % compare_dir


                arg_list = ["compare_files", "%5d" % p, "%5d" % q, "%d" % iframe]
                subprocess.call(arg_list)


def write_ini_files(input_file='create_run.ini'):
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.read(input_file)

    scheduler     = config.get('Defaults','scheduler').partition('#')[0].strip()

    execname      = config.get('Problem','execname').partition('#')[0].strip()

    # This is for estimating total time.  Only time for the uniform calculation is
    # estimated;  adaptive time is assumed to be much less.
    time_per_grid = float(config.get('Problem','time_per_grid').partition('#')[0].strip())

    # Supply a baseline dt and associated grid resolution so we can scale the
    # time step appropriately for different resolutions.
    dt_eff_res    = int(config.get('Problem','dt_eff_res').partition('#')[0].strip())
    dt_fixed      = float(config.get('Problem','dt_fixed').partition('#')[0].strip())

    # Uniform or adaptiev
    mode    = config.get('Run','mode').partition('#')[0].strip()

    # weak or strong
    scaling = config.get('Run','scaling').partition('#')[0].strip()

    # Subcycling (used only in adaptive case)
    sc = config.get('Run','subcycle').partition('#')[0].strip()
    subcycle = sc in ['T','True','1']

    try:
        nw = config.get('Run','noweightedp').partition('#')[0].strip()
        noweightedp = nw in ['T','True','1']
    except:
        noweightedp = False

    # scale-uniform : Scale uniform calculations to save time
    try:
        sw = config.get('Run','scale-uniform').partition('#')[0].strip()
        scale_uniform = sw in ['T','True','1']
    except:
        scale_uniform = False

    mx0     = int(config.get('Run', 'mx').partition('#')[0].strip())
    minlevel0  = int(config.get('Run', 'minlevel').partition('#')[0].strip())
    maxlevel0  = int(config.get('Run', 'maxlevel').partition('#')[0].strip())
    proc0   = int(config.get('Run','proc').partition('#')[0].strip())
    tfinal0 = float(config.get('Run','tfinal').partition('#')[0].strip())

    # nbjobs : determines length of processor sequence, e.g. [1,4,16,64,...]
    njobs   = int(config.get('Run','njobs').partition('#')[0].strip())


    # ----------------------------------------
    # Other inputs needed by the user
    # ----------------------------------------

    # Figure out dt needed for first run in this series
    eff_res0 = mx0*2**minlevel0
    dt0 = dt_fixed/(float(eff_res0)/float(dt_eff_res))

    tol = 1e-15
    nout0 = tfinal0/dt0
    if abs(nout0 - np.round(nout0)) > tol:
        print "nout is not an integer; nout = %12.4f\n" % nout0
        sys.exit()
    nout0 = int(nout0)

    nstep0 = nout0/tfinal0
    if abs(nstep0 - np.round(nstep0)) > tol:
        print "nstep is not an integer; nstep = %12.4f; nout = %12.4f\n" % (nstep0,nout0)
        print "Setting nstep to nout"
        nstep0 = int(nout0)
    else:
        nstep0 = int(nstep0)


    # ----------------------------------------
    # Everything should be automatic from here
    # Four types of runs ;
    #   -- uniform with weak scaling
    #   -- uniform with strong scaling
    #   -- adaptive with weak scaling
    #   -- adaptive with strong scaling
    # ----------------------------------------

    R = np.arange(0,njobs)
    procs = proc0*4**R     # Key idea : Run on different processors


    mx       = np.empty(R.shape)
    nout     = np.empty(R.shape)
    nstep    = np.empty(R.shape)
    minlevel = np.empty(R.shape)
    maxlevel = np.empty(R.shape)
    tfinal   = np.empty(R.shape)
    dt       = np.empty(R.shape)


    # These are fixed for all possible runs.
    mx.fill(mx0)
    nout.fill(nout0)
    nstep.fill(nstep0)

    if mode == 'uniform':
        if scaling == 'weak':
            # Quadruple proc count :
            #    -- Double effective resolution
            #    -- cut dt in half
            # With "scale_uniform" set:
            #    -- cut tfinal in half
            #    -- nout and nstep don't change
            # If "scale_uniform" not set:
            #    -- tfinal is constant
            #    -- nout and nstep are doubled
            minlevel = minlevel0 + R    # Increase minlevel=maxlevel
            dt = dt0/2**R
            if scale_uniform:
                tfinal = tfinal0/2**R       # saves time
            else:
                tfinal.fill(tfinal0)
                nout = nout*2**R
                nstep = nstep*2**R

        elif scaling == 'strong':
            # Quadruple proc count --> nothing changes in problem size
            # with "scale_uniform" set :
            #    -- cut tfinal in half (cut nout/nstep in half as well)
            #    -- scale time later (in plotting, for example)
            minlevel.fill(minlevel0)
            tfinal.fill(tfinal0)   # We will multiply by procs when plotting
            dt.fill(dt0)
            if scale_uniform:
                tfinal = tfinal*procs/np.max(procs)
                nout = nout*procs/np.max(procs)
                nstep = nstep*procs/np.max(procs)

        maxlevel = minlevel   # Always true for uniform case

    elif mode == 'adapt':

        # minlevel is fixed
        minlevel.fill(minlevel0)

        # Run entire simulation, to get true measure of adaptivity
        tfinal.fill(tfinal0)
        dt.fill(dt0)    # Coarse grid dt; is scaled to finer levels.

        if scaling == 'weak':
            # Quadruple processor count
            #   -- double effective resolution
            #   -- dt, nout, nstep are relative to coarsest level
            #      and so remain unchanged.
            maxlevel = minlevel0 + R

        elif scaling == 'strong':
            # Quadruple processor count
            #    -- run same simulation
            maxlevel.fill(maxlevel0)

    nout_uniform = nout*2**(maxlevel-minlevel)
    num_grids = (2**maxlevel)**2  # Total number of uniform grids on all procs
    t = (num_grids*time_per_grid*nout_uniform)/procs  # Assuming ideal scaling
    eff_res = mx*(2**maxlevel)
    grids_per_proc = (2**maxlevel)**2/procs

    # ------------------------------------------
    # Start creating files.
    # ------------------------------------------
    fmt_str_numeric = "# %6d %4d %7d %7d %6d %6d %12.4e %12.4e %8d %12d %8.1f"
    fmt_str_header  = "# %6s %4s %7s %7s %6s %6s %12s %12s %8s %12s %8s"
    tuple_str = ('p','mx','minlev','maxlev','nout','nstep', 'dt','tfinal',
                 'eff_res','grids/proc','t')
    header_line = "# " + "-"*99

    # Output to console
    print header_line
    print fmt_str_header % tuple_str
    print header_line

    # Output to job files (jobs.sh)
    jobfile = open('jobs.sh','w')
    jobfile.write("#!/usr/bin/env bash\n\n")
    jobfile.write(header_line)
    jobfile.write("\n")
    jobfile.write((fmt_str_header + "\n") % tuple_str)
    jobfile.write(header_line)
    jobfile.write("\n")

    cpu_hours = 0
    for i,p in enumerate(procs):
        if grids_per_proc[i] < 4:
            print "Too few grids per proc. No more output files will be created."
            sys.exit()

        prt_tuple = (p,mx[i],minlevel[i],maxlevel[i],nout[i],nstep[i],dt[i],
                     tfinal[i],eff_res[i], grids_per_proc[i],t[i])

        # Print to console and jobs file
        print fmt_str_numeric % prt_tuple
        jobfile.write((fmt_str_numeric + "\n") % prt_tuple)

        if scheduler == 'll':
            cpu_hours = cpu_hours + np.max([p,32])*t[i]/(3600.0)
        else:
            cpu_hours = cpu_hours + t[i]/3600.0

        inifile = "ex_%05d.ini" % p
        ini_file = open(inifile,'w')
        ini_file.write("[user]\n")
        ini_file.write("     example = 0      # no mapping\n")
        ini_file.write("\n")
        ini_file.write("[Options]\n")
        ini_file.write("# Grid dimensions\n")
        ini_file.write("    mx = %d\n" % mx[i])
        ini_file.write("    my = %d\n" % mx[i])
        ini_file.write("\n")
        ini_file.write("    minlevel = %d\n" % minlevel[i])
        ini_file.write("    maxlevel = %d\n" % maxlevel[i])

        if minlevel[i] == maxlevel[i]:
            ini_file.write("    regrid_interval = 0\n")
        else:
            ini_file.write("    regrid_interval = 1\n")
            ini_file.write("    smooth-refine = T\n")
            ini_file.write("    smooth-level = %d\n" % (maxlevel[i]-1))

        if noweightedp:
            ini_file.write("noweightedp = T\n")
        else:
            ini_file.write("noweightedp = F\n")

        ini_file.write("\n")
        ini_file.write("    tfinal = %12.4f\n" % tfinal[i])
        ini_file.write("    use_fixed_dt = T\n")
        ini_file.write("    initial_dt = %16.8e\n" % dt[i])
        ini_file.write("\n")
        ini_file.write("    outstyle = 3\n")
        ini_file.write("    nout = %d\n" % nout[i])
        ini_file.write("    nstep = %d\n" % nstep[i])
        if mode == 'adapt':
            if subcycle:
                ini_file.write("    subcycle = T\n")
            else:
                ini_file.write("    subcycle = F\n")

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
            rpn = np.min([32,p])
            bgsize = np.ceil(p/(32.0*rpn))*32
            proc_file.write("#@ bg_size = %d\n" % bgsize)
            proc_file.write("#@ queue\n")
            proc_file.write("\n")
            proc_file.write(("runjob --ranks-per-node %d --np %d : " \
                             + "/homec/hbn26/hbn263/projects/forestclaw-build-alt/local/bin/torus " \
                             + "--inifile=ex_%05d.ini\n") %
                            (np.min([32,p]),p,p))
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
                print "PBS scheduler : Unexpected number of processors; p = %d\n" % p
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
            pstr = "print \"Starting process %d with jobid %d on $p " + \
                    "processor(s).\" % (po.pid,jobid)\n"
            s = Template(pstr)
            pstr = s.substitute(p=str(p))
            proc_file.write(pstr)
            proc_file.write("po.wait()\n")
            proc_file.close()


    print header_line
    print "Estimated core-h (hours) : %6.2f (%4.2f%%)" % (cpu_hours,100*cpu_hours/350000.0)


    jobfile.write(header_line)
    jobfile.write("\n")
    jobfile.write("# Estimated core-h (hours) : %6.2f (%4.2f%%)" % \
                  (cpu_hours,100*cpu_hours/350000.0))
    jobfile.write("\n\n");

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

    print ""
    print "Configuration files written : "
    for i,p in enumerate(procs):
        inifile = "ex_%05d.ini" % p
        print "    %s" % inifile

    print ""
    print "jobs submitted from jobs.sh"
    for i,p in enumerate(procs):
        if scheduler == 'll':
            procfile = "p_%05d.ll" % p
        elif scheduler == 'pbs':
            procfile = "p_%05d.pbs" % p
        elif scheduler == 'osx':
            procfile = "p_%05d.py" % p
        print "    %s" % procfile


def launch_jobs(N=1):

    import subprocess
    for i in range(0,N):
        po = subprocess.call(['bash','jobs.sh'])

def compile_results(results_dir=None,results_file='results.out',
                    execname=None):

    if results_dir == None:
        results_dir = os.getcwd()

    # Stuff to grep for
    stats_list = [ 'WALLTIME',
                   'ADVANCE',
                   'Statistics for EXCHANGE',
                   'Statistics for REGRID$',
                   'GHOSTCOMM$',
                   'Statistics for CFL',
                   'Statistics for EXTRA4$',
                   'TIME_PER_GRID',
                   'BUILDREGRID',
                   'BUILDGHOST',
                   'TAGGING']

    # Compile data to this file
    resultsfile = open(results_file,'w')
    resultsfile.write("# " + "-"*152)
    resultsfile.write("\n")
    fstr = "# %6s%8s%6s%6s%6s%8s%12s%12s" + "%12s"*11 + "\n"
    resultsfile.write(fstr % ('jobid','p','mx','min','max','nout','tfinal','grid/proc',
                              'Wall','Advance','Exchange','Regrid','Comm','cfl',
                              'extra','time/grid','buildregrid','buildghost','tagging'))

    resultsfile.write("# " + "-"*152)
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

            if found_bad_file:
                print "Skipping file %s (time exceeded)" % f
                continue

            # jobid
            resultsfile.write("  ")   # Make up for "# " added as comments above
            resultsfile.write("%6d" % jobid)

            # proc count
            resultsfile.write("%8d" % pcount)

            # mx
            for i,l in enumerate(lines):
                if re.search("mx",l):
                    l1 = lines[i].split()
                    mx = int(l1[2])
                    resultsfile.write("%6s" % l1[2])
                    break

            # minlevel
            for i,l in enumerate(lines):
                if re.search("minlevel",l):
                    l1 = lines[i].split()
                    minlevel = int(l1[2])
                    resultsfile.write("%6s" % l1[2])
                    break

            # maxlevel
            for i,l in enumerate(lines):
                if re.search("maxlevel",l):
                    l1 = lines[i].split()
                    maxlevel = int(l1[2])
                    resultsfile.write("%6s" % l1[2])
                    break

            # nout
            for i,l in enumerate(lines):
                if re.search("nout",l):
                    l1 = lines[i].split()
                    nout = int(l1[2])
                    resultsfile.write("%8s" % l1[2])
                    break

            # tfinal
            for i,l in enumerate(lines):
                if re.search("tfinal",l):
                    l1 = lines[i].split()
                    resultsfile.write("%12.2e" % float(l1[2]))
                    break

            # ratio of adaptive grids advanced to uniform grids advanced
            # Account for the fact that nout in this case means something
            # different than in the uniform case.
            nout_uniform = nout*2**(maxlevel-minlevel)
            grids_advanced_uniform = nout_uniform*(2**maxlevel)**2/pcount
            a = []
            for i,l in enumerate(lines):
                if re.search('TIME_PER_GRID',l):
                    a = i
                    break

            missing_value = False
            try:
                l2 = lines[a+2].split()  # Data if two lines below "TIME_PER_GRID"
                tpg = float(l2[5])
            except:
                print "TIME_PER_GRID not found in %s" % f
                missing_value = True

            a = []
            for i,l in enumerate(lines):
                if re.search('WALLTIME',l):
                    a = i
                    break

            try:
                l2 = lines[a+2].split()
                wt = float(l2[5])
            except:
                print "WALLTIME not found in %s" % f
                missing_value = True

            if not missing_value:
                grids_advanced_per_proc = wt/tpg
                adapt_ratio = 100*float(grids_advanced_per_proc/grids_advanced_uniform)
                resultsfile.write("%12d" % np.round(grids_advanced_per_proc))
            else:
                resultsfile.write("%12s" % np.nan)

            # Get everything else in list
            for k,w in enumerate(stats_list):
                for i,l in enumerate(lines):
                    if re.search(w,l):
                        a = i
                        break
                try:
                    l2 = lines[a+2].split()
                    resultsfile.write("%12.4e" % float(l2[5]))
                except:
                    print "%s in file %s not found" % (w,f)
                    resultsfile.write("%12f" % (np.nan))

            resultsfile.write('\n')

            run_file.close()

    resultsfile.close()

# Read in an array of results that can be used for both weak
# and strong scaling.
def read_results_files(dir_list, subdir = None, results_in = None,
                       results_file='results.out',execname=None):
    import re

    # Get data from directory names
    dirs = os.listdir(os.getcwd())

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
                    print "File %s not found in %s" % (results_file,f)
                    continue

            s = s[2].partition('_')
            l = int(s[0])
            levels1.append(l)

            p = int(s[2])
            procs1.append(p)

    t = zip(mx1,procs1,levels1)  # All possible combos

    # Create array of all (procs/levels).  Not every combo will have data.
    mx = list(set(mx1))
    procs = list(set([x[1] for x in t]))
    levels = list(set([x[2] for x in t]))

    # Initialize dictionaries for all values
    jobs = dict.fromkeys(mx)
    for m in mx:
        jobs[m] = dict.fromkeys(set(procs))
        for p in procs:
            jobs[m][p] = dict.fromkeys(set(levels))

    # Initialize dictionaries for all values
    for results_dir in dirs:
        rundir = int(results_dir.partition('_')[2])
        if subdir == None:
            d = results_dir
        else:
            d = os.path.join(results_dir,subdir)

        files = os.listdir(d)
        print files
        for f in files:
            if re.match(pattern,f):

                # Load data file to read in (mx,levels) for this file
                if subdir == None:
                    rf = os.path.join(results_dir,results_file)
                else:
                    rf = os.path.join(results_dir,subdir,results_file)

                if not os.path.exists(rf):
                    jobs[m][p][l] = None
                    continue

                data = np.loadtxt(rf)

                try:
                    # Filter to get data only for this (p,l) combination
                    msk = (data[:,1] == p) & (data[:,4] == l)
                except:
                    print "Something went wrong in reading file %s" % fname
                    continue

                job = {}
                job["jobid"]         = np.array(data[msk,0])
                job["procs"]         = np.array(data[msk,1])
                job["mx"]            = np.array(data[msk,2])
                job["minlevel"]      = np.array(data[msk,3])
                job["maxlevel"]      = np.array(data[msk,4])
                job["nout"]          = np.array(data[msk,5])
                job["tfinal"]        = np.array(data[msk,6])
                job["grids_per_proc"]= np.array(data[msk,7])
                job["walltime"]      = np.array(data[msk,8])
                job["advance"]       = np.array(data[msk,9])
                job["exchange"]      = np.array(data[msk,10])
                job["regrid"]        = np.array(data[msk,11])
                job["ghostcomm"]     = np.array(data[msk,12])
                job["cfl"]           = np.array(data[msk,13])
                job["extra"]         = np.array(data[msk,14])
                job["time_per_grid"] = np.array(data[msk,15])
                job["buildregrid"]   = np.array(data[msk,16])
                job["buildghost"]    = np.array(data[msk,17])
                job["tagging"]    = np.array(data[msk,18])

                jobs[m][p][l] = job

    return jobs

def print_jobs(jobs,val2plot):

    mx = sorted(jobs.keys())

    int_list = ['jobid','procs','mx','minlevel','maxlevel','nout', 'grids_per_proc',
                'rundir']

    for m in mx:
        procs = sorted(jobs[m].keys())
        levels = sorted(jobs[m][procs[0]].keys())
        data = np.empty((len(procs),len(levels)))
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

                try:
                    d = job[val2plot]
                    fmt_int = val2plot in int_list
                except:
                    try:
                        d,fmt_int = val2plot(job,mx=m,proc=p,level=l,all=jobs)
                    except:
                        print "print_results : Invalid 'val2plot'"
                        sys.exit()

                try:
                    data[i,j] = np.average(d)
                except:
                    print "job (mx=%d,proc=%d,level=%d) is empty" % (m,p,l)
                    data[i,j] = np.nan


        # Header
        print ""
        print val2plot
        linelen = 8 + 3 + 12*len(levels)
        line = "-"*linelen
        print line
        mxstr = "mx = %d" % m
        header_format = "{mxstr:>8}" + "{sep:>3}" + "{:>12}" * len(levels)
        print header_format.format(mxstr = mxstr, sep="|",*levels)
        print line
        # Format integer and float values separately
        if (val2plot in int_list) or fmt_int:
            row_format = "{:>8d}" + "{sep:>3}" + "{:>12.0f}"*len(levels)
        else:
            row_format = "{:>8d}" + "{sep:>3}" + "{:>12.2e}"*len(levels)
        for p, d in zip(procs, data):
            print row_format.format(p, sep="|", *d)

        print line


def plot_results(jobs,start_point,val2plot='walltime',
                 scaling='strong',scale_uniform=False):

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    markers = {}
    markers["walltime"] = r'$\clubsuit$'
    markers["advance"] = u's'
    markers["exchange"] = u'h'
    markers["regrid"] = u'v'
    markers["ghostcomm"] = u'^'
    markers["time_per_grid"] = u'o'
    markers["cfl"] = u'p'
    markers["extra"] = u'd'

    colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')

    mx = sorted(jobs.keys())
    ph = []
    mb_data = []

    # each start point generates a plot, either down a column (strong), across a
    # diagonal (weak) or from back to front (block size comparison)
    tlist = []
    for i,sp in enumerate(start_point):
        m = sp[0]
        p = sp[1]
        l = sp[2]
        job = jobs[m][p][l]

        # "headers" for this mx table
        procs = sorted(jobs[m].keys())
        levels = sorted(jobs[m][p].keys())

        if scaling  == 'weak':  # Work down/right across diagonal
            mx_sub = [sp[0] for p in procs if p >= sp[1]]
            procs_sub =  [p for p in procs if p >= sp[1]]
            levels_sub = [l for l in levels if l >= sp[2]]

        elif scaling == 'strong':  # Work down colum
            # assume that mx is fixed
            mx_sub = [sp[0] for p in procs if p >= sp[1]]
            procs_sub = [p for p in procs if p >= sp[1]]
            levels_sub = [sp[2] for p in procs if p >= sp[1]]

        elif scaling == 'block':  # Work down/left and down mx
            mx_sub = [m for m in mx if m >= sp[0]]
            levels_sub = [l for l in levels if l <= sp[2]]
            levels_sub.reverse()
            procs_sub = [sp[1] for l in levels]

        elif scaling == 'resolution':
            # scaling as resolution is increased
            levels_sub = [l for l in levels if l >= sp[2]]
            mx_sub = [sp[0] for l in levels_sub]
            procs_sub = [sp[1] for l in levels_sub]

        elif scaling == 'superlinear':
            # Decrease grid size while increasing proc count; fixed resolution
            mx_sub = [m for m in mx if m <= sp[0]]
            mx_sub.reverse()
            # This assumes that jobs[8] has the most complete list of levels.
            levels_sub = sorted([l for l in jobs[8][1].keys() if l >= sp[2]])
            procs_sub = [p for p in procs if p >= sp[1]]

        t = zip(mx_sub,procs_sub,levels_sub)

        # Internal plotting routine
        phandle,mb = plot_results_internal(val2plot,t,jobs,scaling,markers,colors,
                                           scale_uniform)

        phandle.set_color(colors[i])
        plt.draw()
        ph.append(phandle)
        mb_data.append(mb)

        tlist.extend(t)

    # Set up axis
    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%d'))

    if scaling == 'block':
        ax.xaxis.set_major_locator(plt.FixedLocator(mx))
        ax.set_xlabel("Block size",fontsize=16)
        l1 = np.log(np.min(mx))/np.log(2)
        l2 = np.log(np.max(mx))/np.log(2)
        plt.xlim([2**(l1-0.5),2**(l2+0.5) ])

    elif scaling == 'resolution':
        levels = sorted(jobs[8][1].keys())
        ax.xaxis.set_major_locator(plt.FixedLocator(levels))
        ax.set_xlabel("Levels",fontsize=16)
        l1 = np.min(levels)
        l2 = np.max(levels)
        plt.xlim([l1-0.5,l2+0.5])

    else:
        ax.xaxis.set_major_locator(plt.FixedLocator(procs))
        ax.set_xlabel("Processor count",fontsize=16)
        p1 = np.log(np.min(procs))/np.log(4)
        p2 = np.log(np.max(procs))/np.log(4)
        plt.xlim([4**(p1-0.5), 4**(p2+0.5)])
        if scaling == 'weak':
            ax.set_yscale('linear')
            plt.ylim([0,110])
            plt.grid(b=True,which='major')
            plt.grid(b=True,which='minor')
            plt.minorticks_on()

    plt.setp(ax.get_xticklabels(),fontsize=14)
    plt.setp(ax.get_yticklabels(),fontsize=14)

    if scaling == 'weak':
        ax.set_ylabel("%s" % ("grids advanced/time"),fontsize=16)
    else:
        ax.set_ylabel("%s" % ("time"),fontsize=16)

    plt.draw()

    return ph,mb_data


def plot_results_internal(val2plot,t,jobs,scaling,markers,colors,
                          scale_uniform=False):

    import matplotlib.pyplot as plt

    if isinstance(val2plot,list):
        print "Only one val2plot allowed per plot"
        sys.exit()

    v = val2plot

    mode = 'uniform'
    ph = []
    mx = sorted(jobs.keys())
    gpp_avg = []
    y_avg = []
    grids_advanced =[]

    # First collect y values (which may need to be averaged)
    for idx_tuple in t:

        m = idx_tuple[0]
        p = int(idx_tuple[1])
        l = idx_tuple[2]
        job = jobs[m][p][l]
        if job == None:
            y_avg = np.append(y_avg,np.nan)
            continue

        # Each job might include several runs, which must be averaged
        gpp_avg = np.append(gpp_avg,np.round(np.average(job["grids_per_proc"]/job["nout"])))
        grids_advanced = np.append(grids_advanced,
                                   np.round(np.average((job["walltime"]/job["time_per_grid"]))))

        try:
            data = np.average(job[v])
        except:
            try:
                # Call function val2plot
                data = np.average(v(job,mx=m,proc=p,level=l,all=jobs))
            except:
                import pdb
                pdb.set_trace()

                print "Invalid val2plot"
                sys.exit()


        if scaling == 'weak':
            y = grids_advanced[-1]/data
        else:
            y = data

        y_avg = np.append(y_avg,y)

    if scaling == 'weak':
        # Compute efficiency
        y_avg = 100*y_avg/y_avg[0]


    # Plot ideal scaling
    tp = zip([x[0] for x in t],[x[1] for x in t],[x[2] for x in t],y_avg)
    mx     = [x[0] for x in tp if not np.isnan(x[3])]
    procs  = [x[1] for x in tp if not np.isnan(x[3])]
    levels = [x[2] for x in tp if not np.isnan(x[3])]
    y_avg  = [x[3] for x in tp if not np.isnan(x[3])]

    if scaling == 'weak':
        # Ideal scaling is constant
        ideal = [y_avg[0] for x in tp if not np.isnan(x[3])]  # Constant scaling
        plt.loglog(procs,ideal,'k.--',markersize=15)

    elif scaling == 'strong':
        # Ideal scaling should start at p = 16, not p=1
        ideal = y_avg[0]*np.min(procs)/np.array(procs)
        plt.loglog(procs,ideal,'k.--',markersize=15)

    elif scaling == 'superlinear':
        ideal = y_avg[0]*np.min(procs)/np.array(procs)
        plt.loglog(procs,ideal,'k.--',markersize=15)

    elif scaling == 'resolution':
        if v == 'time_per_grid':
            R = np.arange(0,len(y_avg))
            ideal = [y_avg[0] for x in t]
        else:
            R = np.array([i for i,x in enumerate(tp) if not np.isnan(x[3])])
            ideal = y_avg[0]*8**R

        plt.semilogy(levels,ideal,'k.--',markersize=15)

    elif scaling == 'block':
        if v == 'time_per_grid':
            R = np.arange(0,len(y_avg))
            ideal = y_avg[0]*4**R
            plt.loglog(mx,ideal,'k.--',markersize=15)
        else:
            # Try a parabolic fit
            if len(mx) > 2:
                pc = np.polyfit(np.log(mx),np.log(y_avg),2)
                m = np.logspace(np.log10(4),np.log10(256),40)
                plt.loglog(m,np.exp(np.polyval(pc,np.log(m))),'k--')

    # get best-fit slope
    mb = None
    if scaling in ['strong','weak','superlinear']:
        mb = np.polyfit(np.log(procs),np.log(y_avg),1)

    elif scaling == 'resolution':
        nsub = np.array(levels) - levels[0] + 1
        pc = np.polyfit(nsub,np.log(y_avg),1)
        mb = np.array(np.exp(pc))

    elif scaling == 'block':
        if v == 'time_per_grid':
            mx = np.array([x[0] for x in t])
            nsub = [i for i,te in enumerate(t)]
            pc = np.polyfit(nsub,np.log(y_avg),1)
            mb = np.array(np.exp(pc))

    try:
        mv = markers[v]
    except:
        mv = u'o'


    if scaling == 'block':
        ph = plt.loglog(mx,y_avg,marker=mv,markersize=10)
    elif scaling in ['strong','superlinear']:
        ph = plt.loglog(procs,y_avg,marker=mv,markersize=10)
    elif scaling == 'weak':
        ph = plt.loglog(procs,y_avg,marker=mv,markersize=10)
    elif scaling in 'resolution':
        ph = plt.semilogy(levels,y_avg,marker=mv,markersize=10)

    return ph[0],mb


if __name__ == "__main__":
    import sys
    write_ini_files()
