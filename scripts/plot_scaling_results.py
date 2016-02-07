#!/usr/bin/env python

import parallel_tools
import matplotlib.pyplot as plt
import numpy as np
import sys

# --------------------------------------------------------------------------
# This is a sample script that uses the function 'plot_results' in the python
# module 'parallel_tools'.  To use this module, you must have a few
# things in place :
#
#   1.  The current working directory must contain directories with names
#       MMM_LL_PPPPP, where 'MMM' indicates the block size, 'LL' indicates
#       the level and 'PPPPP' indicates the processor count.  Each number
#       is padded with zeros.  Example : '008_04_00016'
#   2.  Directories must contain a results file (default : results.out)
#       created from calling the 'compile_results' function in the parallel_tools
#       module.  This file contains results from jobs run on a sequence of
#       processors, for problems at a fixed effective resolution and block
#       size.
#
# See 'parallel_tools' for more info.
# --------------------------------------------------------------------------

def val2plot(job,mx=None,proc=None,level=None,all=None):

    # v = job["walltime"]/job["advance_steps"]
    v = job["walltime"]

    fmt_int = False

    # Efficiency
    if not all==None and False:
        # Compute the efficiency
        wall = {}
        wall[1] = np.average(jobs[mx][1][level]["walltime"])
        wall[proc] = np.average(jobs[mx][proc][level]["walltime"])
        v = 100*(wall[1]/wall[proc])/float(proc)
        fmt_int = True

    return v, fmt_int


execname = 'Torus'
platform = 'Juqueen'

# Plot data from directories run_XXX, where XXX is a zero-padded digit
# in dir_list.
dir_list = range(6,11)

# ---------------------------------------------------
# read in data files stored in results files
# ---------------------------------------------------

#jobs = parallel_tools.read_results_files(dir_list,results_in = 'results.in')
#parallel_tools.print_jobs(jobs,val2plot)
#sys.exit()

jobs = parallel_tools.read_results_files(dir_list)

# parallel_tools.print_jobs(jobs,val2plot)


# ---------------------------------------------------------------------
# Create plots.
# The 'scaling' parameter determines how the data will be traversed
#
# scaling = 'strong': Fixes the resolution, block size and level.
# and level) and traverses down columns of data.
#
# scaling = 'weak': Fixes the block size, but traverse diagonally down
# across rows and columns simultaneously.
#
# scaling = 'superlinear' : Like strong scaling, with fixed resolution, but
# block size descreases with increase in processor count, fixing the number
# of blocks per processor.  Improved cache performance could lead to
# superlinear scaling.
#
# scaling = 'resolution' : Fixed blocksize and processor count, but increase
# resolution.  Should expect to see 8-fold increase in wall time.
#
# scaling = 'block' : Fix resolution and processor count, but increase block
# size.  One might expect to see slightly parabolic profile as performance of
# cache competes with communication costs.
#
# ---------------------------------------------------------------------------

name = 'Walltime'
v = 'walltime'   # Quantity to plot
scaling = 'strong'

if scaling == 'weak':
    mx = [8]*5
    procs =  [1,4,16]
    levels = [5]*3
    start_points = zip(mx,procs,levels)
    # start_points.reverse()


elif scaling == 'strong':
    # Fix mx and level; vary proc count
    mx = [8]*5
    levels = [5,6,7,8,9]
    procs = [1,1,1,4,16]
    start_points = zip(mx,procs,levels)

elif scaling == 'superlinear':
    start_points = [(128,1,4),(64,1,4),(32,1,4),(16,1,4)]

elif scaling == 'resolution':
    mx = [16]*6
    levels = [3,3,4,5,6,7]
    procs = [1,4,16,64,256,1024]
    start_points = zip(mx,procs,levels)
    # start_points.reverse()

elif scaling == 'block':
    # Block scaling
    mx = [8]*4
    procs = [1]*4
    levels = [5,6,7,8]
    start_points = zip(mx,procs,levels)
    start_points.reverse()

# Plot results and get best fit lines (if they are available)
parallel_tools.print_jobs(jobs,v)
handles,mb,t = parallel_tools.plot_results(jobs,start_points,val2plot=v,
                                         scaling=scaling)

ax = plt.gca()
# plt.yscale('linear')
if scaling == 'weak':
    plt.ylim([-20, 110])



# Set labels for legends and add the legend
for i,phandle in enumerate(handles):
    sp = start_points[i]
    if scaling == 'block':
        if v == 'time_per_grid':
            phandle.set_label('(%d,%d,%d); %d (%6.2f)' %
                              (sp[0],sp[1],sp[2],sp[0]*2**sp[2],mb[i][0]))
        else:
            phandle.set_label('(%d,%d,%d); %d' %
                              (sp[0],sp[1],sp[2],sp[0]*2**sp[2]))

    elif scaling == 'strong':
        job = jobs[sp[0]][sp[1]][sp[2]]
        res = job["mx"]*2**job["maxlevel"]
        phandle.set_label('(%d); (%6.2f)' % (res,mb[i][0]))

    elif scaling == 'weak':
        job = jobs[sp[0]][sp[1]][sp[2]]
        gpp = job["advance_steps"]
        phandle.set_label('(%d grids/proc); (%6.2f)' %
                          (np.round(gpp),mb[i][0]))

    else:
        phandle.set_label('(%d,%d,%d) (%6.2f)' %
                        (sp[0],sp[1],sp[2],mb[i][0]))

plt.legend(bbox_to_anchor=(0, 0, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)

# plt.legend(handles=handles,loc=pos)

# add a title

if scaling == 'block':
    plt.title("%s (%s) : %s (%s scaling on %d processor(s))" % \
              (execname.capitalize(),platform.capitalize(),
               name.capitalize(),scaling,start_points[0][1]))

elif scaling == 'strong':
    plt.title("%s (%s) : %s (%s scaling; mx = %d)" % \
              (execname.capitalize(),platform.capitalize(),
               name.capitalize(),scaling,start_points[0][0]))

elif scaling == 'weak':
    plt.title("%s (%s) : %s (%s scaling; mx = %d)" % \
              (execname.capitalize(),platform.capitalize(),
               name.capitalize(),scaling,start_points[0][0]))
elif scaling == 'resolution':
    plt.title("%s (%s) : %s (%s scaling; mx = %d)" % \
              (execname.capitalize(),platform.capitalize(),
               name.capitalize(),scaling,start_points[0][0]))
else:
    plt.title("%s (%s) : %s (%s scaling)" % \
              (execname.capitalize(),platform.capitalize(),
               name.capitalize(),scaling))

plt.savefig("walltime_weak.png")
plt.show()
