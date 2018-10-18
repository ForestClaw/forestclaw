#!/bin/env python3

import sys
import os
import subprocess
import random

# Redhawk has only 4 cores per node on compute nodes;  2 GPUs per node; one PCI bus

if sys.argv.__len__() == 1:
    print("Usage : ")
    print("   run_radial.py <nodes> <ntasks-per-node> <gpus-per-node> <nout>")
    print("")
    print("          nodes           : Number of nodes to use")
    print("          ntasks-per-node : MPI processes per node")
    print("          gpus-per-node   : GPUs per node (optional)")
    print("          nout            : (optional) Number of output steps")
    print("")
    sys.exit()

nodes          = int(sys.argv[1])
tasks_per_node = int(sys.argv[2])
gpus           = int(sys.argv[3])

cpus_per_task  = 1  # CPU cores per MPI process
mpi_tasks = nodes*tasks_per_node

if (gpus == 0):
    gpustr = "--user:cuda=F"
else:
    gpustr = "--user:cuda=T"

arg_list = ["srun", "--gres=gpu:2".format(gpus),
            "--nodes={:d}".format(nodes),
            "--ntasks-per-node={:d}".format(tasks_per_node),
            "--cpus-per-task={:d}".format(cpus_per_task),            
            "radial", gpustr]     

np = mpi_tasks
jobid = random.randint(1000,9999)
outfile = "radial_{:05d}.o{:d}".format(np,jobid)
f = open(outfile,'w')
po = subprocess.Popen(arg_list,stdout=f)
print("Starting process {:d} with jobid {:d} on {:d} processor(s).".format(po.pid,jobid,np))
#po.wait()
