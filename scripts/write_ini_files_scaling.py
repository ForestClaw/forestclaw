import os
import sys
import numpy as np
from string import Template
import ConfigParser

# ----------------------------------------------------------
# Create a series of .ini files for determining scaling.
#    -- Assume that each file is for a different processor
#    -- files names created use proc. number
# ----------------------------------------------------------


# ---------------------
# Basic inputs
# ---------------------

config = ConfigParser.SafeConfigParser(allow_no_value=True)
config.read('create_run.ini')

scheduler     = config.get('Defaults','scheduler').partition('#')[0].strip()

execname      = config.get('Problem','execname').partition('#')[0].strip()

# This is for estimating total time.  Only time for the uniform calculation is
# estimated;  adative time is assumed to be much less.
time_per_grid = float(config.get('Problem','time_per_grid').partition('#')[0].strip())

# Supply a baseline dt and associated grid resolution so we can scale the
# time step appropriately for different resolutions.
dt_eff_res    = int(config.get('Problem','dt_eff_res').partition('#')[0].strip())
dt_fixed      = float(config.get('Problem','dt_fixed').partition('#')[0].strip())

mode    = config.get('Run','mode').partition('#')[0].strip()
scaling = config.get('Run','scaling').partition('#')[0].strip()
subcycle = bool(config.get('Run','subcycle').partition('#')[0].strip())

mx0     = int(config.get('Run', 'mx').partition('#')[0].strip())
minlevel0  = int(config.get('Run', 'minlevel').partition('#')[0].strip())
maxlevel0  = int(config.get('Run', 'maxlevel').partition('#')[0].strip())
proc0   = int(config.get('Run','proc').partition('#')[0].strip())
njobs   = int(config.get('Run','njobs').partition('#')[0].strip())

tfinal0 = float(config.get('Run','tfinal').partition('#')[0].strip())


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
# ----------------------------------------
R = np.array(range(0,njobs))
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
        # Double effective resolution
        minlevel = minlevel0 + R    # minlevel == maxlevel
        tfinal = tfinal0/2**R   # saves time
        dt = dt0/2**R

    elif scaling == 'strong':
        # Solve the same problem on increasing number of procs.
        minlevel.fill(minlevel0)
        tfinal.fill(tfinal0)
        dt.fill(dt0)

    maxlevel = minlevel   # Always true for uniform case

elif mode == 'adapt':

    # minlevel is fixed
    minlevel.fill(minlevel0)

    # Run entire simulation, to get true measure of adaptivity
    tfinal.fill(tfinal0)
    dt.fill(dt0)    # Coarse grid dt; is scaled to finer levels.

    if scaling == 'weak':
        # Double effective resolution
        maxlevel = minlevel0 + R

    elif scaling == 'strong':
        # Solve same problem on increasing number of procs
        maxlevel.fill(maxlevel0)

# ------------------------------------------
# Start creating files.
# ------------------------------------------
cpu_hours = 0
fmt_str_numeric = "# %6d %4d %7d %7d %6d %6d %12.4e %12.4e %8d %12d %8.1f"
fmt_str_header  = "# %6s %4s %7s %7s %6s %6s %12s %12s %8s %12s %8s"
tuple_str = ('p','mx','minlev','maxlev','nout','nstep', 'dt','tfinal','eff_res','grids/proc','t')

# Output to console
print "# " + "-"*99
print fmt_str_header % tuple_str
print "# " + "-"*99

# Output to job files (jobs.sh)
jobfile = open('jobs.sh','w')
jobfile.write("#!/usr/bin/env bash\n\n")
jobfile.write("# " + "-"*99)
jobfile.write("\n")
jobfile.write((fmt_str_header + "\n") % tuple_str)
jobfile.write("# " + "-"*99)
jobfile.write("\n")

for i,p in enumerate(procs):
    num_grids = (2**maxlevel[i])**2  # Number of uniform grids
    t = (num_grids*time_per_grid*nout[i])/p  # Assuming ideal scaling
    eff_res = mx[i]*(2**maxlevel[i])
    grids_per_proc = (2**maxlevel[i])**2/p
    if grids_per_proc < 4:
        print "Too few grids per proc. No more output files will be created."
        sys.exit()

    prt_tuple = (p,mx[i],minlevel[i],maxlevel[i],nout[i],nstep[i],dt[i],
                 nout[i]*dt[i],eff_res, grids_per_proc,t)

    # Print to console and jobs file
    print fmt_str_numeric % prt_tuple
    jobfile.write((fmt_str_numeric + "\n") % prt_tuple)

    if scheduler == 'll':
        cpu_hours = cpu_hours + np.max([p,32])*t/(3600.0)
    else:
        cpu_hours = cpu_hours + t/3600.0

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
        proc_file.write("runjob --ranks-per-node %d --np %d : %s --inifile=ex_%05d.ini\n" %
                        (np.min([32,p]),p,execname,p))
        proc_file.close()

    elif scheduler == 'pbs':
        # Kestrel (Linux Cluster)
        procfile = "p_%05d.pbs" % p
        proc_file = open(procfile,'w')
        proc_file.write("#!/bin/bash\n")
        proc_file.write("#PBS -l walltime=00:30:00\n")
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
        proc_file.write("mpirun $PBS_O_WORKDIR/torus --inifile=%s\n" % inifile)

    elif scheduler == 'osx':
        # Mac (laptop or desktop). Write a Python script that calls
        # mpirun as a subprocess.
        procfile = "p_%05d.py" % p
        proc_file = open(procfile,'w')
        proc_file.write("# comment = \"%s example : eff. resolution = %d x %d\"\n" % \
                        (execname.capitalize(),eff_res,eff_res))

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
        pstr = "print \"Starting process %d with jobid %d on $p processor(s).\" % (po.pid,jobid)\n"
        s = Template(pstr)
        pstr = s.substitute(p=str(p))
        proc_file.write(pstr)
        proc_file.write("po.wait()\n")
        proc_file.close()


print "# " + "-"*99
print "Estimated core-h (hours) : %6.2f (%4.2f%%)" % (cpu_hours,100*cpu_hours/350000.0)


jobfile.write("# " + "-"*99)
jobfile.write("\n")
jobfile.write("# Estimated core-h (hours) : %6.2f (%4.2f%%)" % (cpu_hours,100*cpu_hours/350000.0))
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
