import sys
import os
import subprocess
import random

nprocs = 1
arg_list = ["mpirun","-n",str(nprocs),"replicated","--inifile", \
    "ex_{:05d}.ini".format(nprocs)]
jobid = random.randint(1000,9999)
outfile = "replicated_{:05d}.o{:d}".format(nprocs,jobid)
f = open(outfile,'w')
po = subprocess.Popen(arg_list,stdout=f)
print("Starting process %d with jobid %d on %d processor(s)." % (po.pid,
                                                                 nprocs,jobid))
# po.wait()
