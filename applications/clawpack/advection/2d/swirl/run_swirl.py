import sys
import os
import subprocess
import random

# srun -o swirl_00004.o%j --gres=gpu:1 --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 swirl_mpi --

# run_swirl.py <nodes> <gpus-per-node> <ntasks-per-node> <cpus-per-task> <use_cuda>

if sys.argv.__len__() != 5:
    print("Usage : ")
    print("   run_swirl.py <nodes> <gpus-per-node> <ntasks-per-node> <use_cuda>")
    print("")
    print("          nodes           : Number of nodes to use")
    print("          gpus-per-node   : GPUs per node")
    print("          ntasks-per-node : MPI processes per node")
    print("")
    sys.exit()

nodes          = int(sys.argv[1])
gpus           = int(sys.argv[2])
tasks_per_node = int(sys.argv[3])
use_cuda       = sys.argv[4]
cpus_per_task  = 1  # CPU cores per MPI process

mpi_tasks = nodes*tasks_per_node

arg_list = ["srun", 
            "--gres=gpu:{:d}".format(gpus), 
            "--nodes={:d}".format(nodes),
            "--ntasks-per-node={:d}".format(tasks_per_node),
            "--cpus-per-task={:d}".format(cpus_per_task),            
            "swirl", "--user:cuda={:s}".format(use_cuda)]

np = mpi_tasks
jobid = random.randint(1000,9999)
outfile = "swirl_0000{:d}.o{:d}".format(np,jobid)
f = open(outfile,'w')
po = subprocess.Popen(arg_list,stdout=f)
print("Starting process {:d} with jobid {:d} on {:d} processor(s).".format(po.pid,jobid,np))
#po.wait()
