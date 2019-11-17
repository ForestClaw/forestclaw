# comment = "Multigrid example"

import sys
import os
import subprocess
import random


if sys.argv.__len__() == 1:
    print("Usage")
    print("      python p_00001.py  <mx> <minlevel> <maxlevel>")
    print("")
elif sys.argv.__len__() == 4:
    mx = my = int(sys.argv[1])
    minlevel = int(sys.argv[2])
    maxlevel = int(sys.argv[3])
else:
    print("Please enter maxlevel")
    sys.exit()

np = 1
exec = "mgtest"

arg_list = ["mpirun","-n",str(np),exec,
    "--clawpatch:mx={:d}".format(mx),
    "--clawpatch:my={:d}".format(my),
    "--minlevel={:d}".format(minlevel),
    "--maxlevel={:d}".format(maxlevel),
    "--report-timing-verbosity=all"]

jobid = random.randint(1000,9999)
outfile = "{:s}_0000{:d}.o{:d}".format(exec,np,jobid)
f = open(outfile,'w')
po = subprocess.Popen(arg_list,stdout=f)
print("Starting process {:d} with jobid {:d} on {:d} processor(s).".format(po.pid,jobid,np))
#po.wait()
