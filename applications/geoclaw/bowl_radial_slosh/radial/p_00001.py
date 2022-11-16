# comment = "Torus example : eff. resolution = 2048 x 2048"
import sys
import os
import subprocess
import random

arg_list = ["mpirun","-n","1","bowl", "--verbosity=production"]
jobid = random.randint(1000,9999)
outfile = "bowl_00001.o{:d}".format(jobid)
f = open(outfile,'w')
po = subprocess.Popen(arg_list,stdout=f)
print("Starting process {:d} with jobid {:d} on 1 processor(s).".format(po.pid,jobid))
po.wait()
