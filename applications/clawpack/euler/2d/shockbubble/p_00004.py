# comment = "Shockbubble example : eff. resolution = 2048 x 512"
import sys
import os
import subprocess
import random

np = 4
exec = "shockbubble"
options_file = 'fclaw_options.ini'

arg_list = ["mpirun","-n",str(np),exec,"-F",options_file]
jobid = random.randint(1000,9999)
outfile = "{:s}_0000{:d}.o{:d}".format(exec,np,jobid)
f = open(outfile,'w')
po = subprocess.Popen(arg_list,stdout=f)
print("Starting process {:d} with jobid {:d} on {:d} processor(s).".format(po.pid,jobid,np))
# po.wait()
