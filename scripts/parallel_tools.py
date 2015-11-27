import os
import sys
import numpy as np
from string import Template
import re
import ConfigParser

#  --------------------------
#  Four scaling categories :
#  --------------------------
#
#  Uniform (weak) (minlevel = maxlevel)
#  levels/procs  |  5    6    7    8
#  --------------|-------------------
#    1  proc     |  W
#    4  procs    |       W
#    16 procs    |            W
#    64 procs    |                W    (possibly scale to save time?)
#
#
#  Uniform (strong) (minlevel = maxlevel)
#  levels/procs  |  5    6    7    8
#  --------------|-------------------
#    1  proc     |                 S   (possibly scale to save time?)
#    4  procs    |                 S
#    16 procs    |                 S
#    64 procs    |                 S
#
#
#  (adaptive, weak)  (minlevel <= maxlevel)
#  levels/procs  |  5    6    7    8
#  --------------|-------------------
#    1  proc     |  W
#    4  procs    |  W    W
#    16 procs    |  W    W    W
#    64 procs    |  W    W    W    W
#
#
#  (adaptive, strong)  (midlevel < maxlevel)
#  levels/procs  |  5    6    7    8
#  --------------|-------------------
#    1  proc     |  S    S    S    S
#    4  procs    |  S    S    S    S
#    16 procs    |  S    S    S    S
#    64 procs    |  S    S    S    S




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
        sys.exit()
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
                            (execname.capitalize(),eff_res,eff_res))
            proc_file.write("#@ error = $(job_name).$(jobid).out\n")
            proc_file.write("#@ output = $(job_name).$(jobid).out\n")
            proc_file.write("\n")
            proc_file.write("#@ environment = COPY_ALL\n")
            if p <= 1024:
                proc_file.write("#@ wall_clock_limit = 00:30:00\n")
            else:
                proc_file.write("#@ wall_clock_limit = 01:00:00\n")

            proc_file.write("#@ notification = error\n")
            proc_file.write("#@ notify_user = donnacalhoun@boisestate.edu\n")
            proc_file.write("#@ job_type = bluegene\n")
            rpn = np.min([32,p])
            bgsize = np.ceil(p/(32.0*rpn))*32
            proc_file.write("#@ bg_size = %d\n" % bgsize)
            proc_file.write("#@ queue\n")
            proc_file.write("\n")
            proc_file.write(("runjob --ranks-per-node %d --np %d : " + \
                             "%s --inifile=ex_%05d.ini\n") %
                            (np.min([32,p]),p,execname,p))
            proc_file.close()

        elif scheduler == 'pbs':
            # Kestrel (Linux Cluster)
            procfile = "p_%05d.pbs" % p
            proc_file = open(procfile,'w')
            proc_file.write("#!/bin/bash\n")
            trun = 2*t[i]
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
        procfile = "p_%05d.py" % p
        print "    %s" % procfile


def launch_jobs(N=1):

    import subprocess
    for i in range(0,N):
        po = subprocess.call(['bash','jobs.sh'])

def compile_results(config_file='create_run.ini'):
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.read(config_file)

    execname  = config.get('Problem','execname').partition('#')[0].strip()
    njobs     = int(config.get('Run','njobs').partition('#')[0].strip())
    proc0     = int(config.get('Run','proc').partition('#')[0].strip())
    mode      = config.get('Run','mode').partition('#')[0].strip()
    sc        = config.get('Run','subcycle').partition('#')[0].strip()
    subcycle  = sc in ['True','T','1']

    R = np.array(range(0,njobs))
    procs = proc0*4**R

    # ----------------------------------------------------------------------------------------
    stats_list = [ 'WALLTIME', 'ADVANCE', 'Statistics for EXCHANGE','Statistics for REGRID$',
                   'GHOSTCOMM$','Statistics for CFL','Statistics for EXTRA4$','TIME_PER_GRID']
    # ----------------------------------------------------------------------------------------

    output_files = os.listdir(os.getcwd())

    resultsfile = open('results.out','w')
    resultsfile.write("# " + "-"*138)
    resultsfile.write("\n")
    fstr = "# %8s%8s%6s%8s%16s" + "%12s"*8 + "\n"
    resultsfile.write(fstr % ('jobid','p','mx','nout','ratio (%)','Wall','Advance',
                              'Exchange','Regrid','Comm','cfl','extra','time/grid'))

    resultsfile.write("# " + "-"*138)
    resultsfile.write("\n")
    for j,p in enumerate(procs):
        fname = '%s_%05d' % (execname,p)
        pcount = procs[j]

        for f in output_files:
            if (f.startswith(fname)):

                # Extract proc count and jobid from output file name
                # E.g. torus_01.o4578
                s = f.split('_')[1].split('.o')
                jobid = int(s[1])

                run_file = open(f,'r')
                lines = run_file.readlines()

                # jobid
                resultsfile.write("  ")   # Make up for "# " added as comments above
                resultsfile.write("%8d" % jobid)

                # proc count
                resultsfile.write("%8d" % pcount)

                # mx
                for i,l in enumerate(lines):
                    if re.search("mx",l):
                        l1 = lines[i].split()
                        mx = int(l1[2])
                        resultsfile.write("%6s" % l1[2])
                        break

                # nout
                for i,l in enumerate(lines):
                    if re.search("nout",l):
                        l1 = lines[i].split()
                        nout = int(l1[2])
                        resultsfile.write("%8s" % l1[2])
                        break

                # minlevel
                for i,l in enumerate(lines):
                    if re.search("minlevel",l):
                        l1 = lines[i].split()
                        minlevel = int(l1[2])
                        break

                # maxlevel
                for i,l in enumerate(lines):
                    if re.search("maxlevel",l):
                        l1 = lines[i].split()
                        maxlevel = int(l1[2])
                        break

                # ratio of adaptive grids advanced to uniform grids advanced
                # Account for the fact that nout in this case means something
                # different than in the uniform case.
                nout_uniform = nout*2**(maxlevel-minlevel)
                grids_advanced_uniform = nout_uniform*(2**maxlevel)**2/pcount
                for i,l in enumerate(lines):
                    if re.search('TIME_PER_GRID',l):
                        a = i
                        break

                l2 = lines[a+2].split()
                tpg = float(l2[5]);

                for i,l in enumerate(lines):
                    if re.search('WALLTIME',l):
                        a = i
                        break
                l2 = lines[a+2].split()
                wt = float(l2[5])

                grids_advanced_per_proc = wt/tpg
                adapt_ratio = 100*float(grids_advanced_per_proc/grids_advanced_uniform)
                resultsfile.write("%16d" % int(grids_advanced_per_proc))


                # Get everything else in list
                for k,w in enumerate(stats_list):
                    for i,l in enumerate(lines):
                        if re.search(w,l):
                            a = i
                            break

                    l2 = lines[a+2].split()
                    resultsfile.write("%12.4e" % float(l2[5]))

                resultsfile.write('\n')


                run_file.close()

    resultsfile.close()

def plot_results(val2plot = 'walltime',config_files='create_run.ini',
                 figname='scaling.png',results_files='results.out',
                 markers=None,colors=None):

    import matplotlib.pyplot as plt

    if not isinstance(results_files,list):
        results_files = [results_files]

    if not isinstance(config_files,list):
        config_files = [config_files]

    jobs = []
    for i,r in enumerate(results_files):

        config = ConfigParser.SafeConfigParser(allow_no_value=True)
        config.read(config_files[i])

        execname  = config.get('Problem','execname').partition('#')[0].strip()
        mode      = config.get('Run','mode').partition('#')[0].strip()
        scaling   = config.get('Run','scaling').partition('#')[0].strip()
        mx = int(config.get('Run','mx').partition('#')[0].strip())
        maxlevel = int(config.get('Run','maxlevel').partition('#')[0].strip())

        try:
            sw = config.get('Run','scale-uniform').partition('#')[0].strip()
            scale_uniform = sw in ['T','True','1']
        except:
            scale_uniform = False

        data = np.loadtxt(r)

        job = {}
        job["jobid"]         = np.array(data[:,0])
        job["procs"]         = np.array(data[:,1])
        job["mx"]            = np.array(data[:,2])
        job["nout"]          = np.array(data[:,3])
        job["grids"]         = np.array(data[:,4])
        job["walltime"]      = np.array(data[:,5])
        job["advance"]       = np.array(data[:,6])
        job["exchange"]      = np.array(data[:,7])
        job["regrid"]        = np.array(data[:,8])
        job["ghostcomm"]     = np.array(data[:,9])
        job["cfl"]           = np.array(data[:,10])
        job["extra"]         = np.array(data[:,11])
        job["time_per_grid"] = np.array(data[:,12])

        job["execname"] = execname
        job["maxlevel"] = maxlevel
        job["mode"] = mode
        job["scaling"] = scaling
        job["mx"] = mx
        job["scale-uniform"] = scale_uniform

        jobs.append(job)

    plot_results_internal(val2plot,jobs,markers,colors)




def plot_results_internal(val2plot,jobs,markers,colors):

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    if not isinstance(val2plot,list):
        val2plot = [val2plot]

    if markers == None:
        markers = []
        for m in Line2D.markers:
            try:
                if len(m) == 1 and m != ' ':
                    markers.append(m)
            except TypeError:
                pass

        styles = markers + [
            r'$\lambda$',
            r'$\bowtie$',
            r'$\circlearrowleft$',
            r'$\clubsuit$',
            r'$\checkmark$']

    if colors == None:
        colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')


    ph = []
    for j,job in enumerate(jobs):

        procs = job["procs"]
        procs_unique = np.unique(job["procs"])

        njobs = job["procs"].shape[0]
        mode = job["mode"]
        scaling = job["scaling"]
        mx = job["mx"]
        maxlevel = job["maxlevel"]
        scale_uniform = job["scale-uniform"]

        grids_advanced = np.round(job["walltime"]/job["time_per_grid"])

        pmax = np.max(procs_unique)
        for i,v in enumerate(val2plot):
            R = np.arange(0,len(procs))
            if scaling == 'weak':
                y_all = job[v]/grids_advanced
            elif scaling == 'strong':
                if mode == 'uniform' and scale_uniform:
                    y_all = job[v]*pmax/procs
                else:
                    y_all = job[v]

            y_avg = []
            y_err = []
            print "%s (%s)" % (v,mode)
            for p in procs_unique:
                y_unique = y_all[np.where(procs == p)]
                y_avg.append(np.average(y_unique))
                y_err.append(np.std(y_unique))
                print "    %12.4e %12.4e" % (y_avg[-1], y_err[-1])
            print ""

            # Get ideal scaling for this case
            p1 = int(np.log(procs_unique[0])/np.log(4.0))
            R = np.arange(p1,len(procs_unique)+p1)
            if scaling == 'weak':
                ideal = np.empty(len(y_avg),)
                ideal.fill(y_avg[0])

            elif scaling == 'strong':
                if v == 'ghostcomm':
                    idx = np.where(procs_unique > 1)[0]
                    ideal = y_avg[idx[0]]/(procs_unique/procs_unique[idx[0]])
                else:
                    ideal = y_avg[0]/(procs_unique/procs_unique[0])

            plt.errorbar(procs_unique,y_avg,yerr = y_err,fmt='o',hold=True)
            ax = plt.gca();
            ax.set_yscale('log')
            phandle = plt.loglog(procs_unique,y_avg,marker=markers[j],
                                 color=colors[j],markersize=10,hold=True,label=v)
            ph.append(phandle[0])

            mb = np.polyfit(np.log(procs_unique),np.log(y_avg),1)
            ph[j].set_label((v + " (%s, slope = %4.2f)") % (mode,mb[0]))

            plt.loglog(procs_unique,ideal,'k.--',markersize=15)


    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%d'))
    ax.xaxis.set_major_locator(plt.FixedLocator(procs_unique))
    plt.setp(ax.get_xticklabels(),fontsize=14)
    plt.setp(ax.get_yticklabels(),fontsize=14)
    ax.set_xlabel("Processor count",fontsize=16)
    if scaling == 'weak':
        plt.legend(handles=ph,labels=val2plot,loc=2)
    else:
        plt.legend(handles=ph,labels=val2plot,loc=3)
    if mode == 'weak':
        ax.set_ylabel("%s" % ("time/grids advanced"),fontsize=16)
    if mode == 'strong':
        ax.set_ylabel("%s" % ("time"),fontsize=16)

    plt.xlim([4**(-0.5), 4**(3.4)])

    # add a title
    eff_res = mx*(2**maxlevel)
    if (len(val2plot) == 1):
        plt.title("%s : %d x %d (%d x %d, %s)" % \
                  (val2plot[0].capitalize(),eff_res,eff_res,mx,mx,mode))
    else:
        plt.title("%s scaling results : %d x %d (%d x %d,%s)" % \
                  (scaling.capitalize(),eff_res,eff_res,mx,mx,mode))

    plt.savefig("scaling.png")
    plt.show()


if __name__ == "__main__":
    import sys
    write_ini_files()
