# comment = "Slotted_disk example : eff. resolution = 1536 x 1536"
import sys
import os
import subprocess
import random

arg_list = ["mpirun","-n","4","slotted_disk","--inifile=sisc.ini"]
jobid = random.randint(1000,9999)
outfile = "slotted_disk_00004.o%d" % (jobid)
f = open(outfile,'w')
po = subprocess.Popen(arg_list,stdout=f)
print "Starting process %d with jobid %d on 4 processor(s)." % (po.pid,jobid)
po.wait()
