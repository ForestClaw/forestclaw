#!/bin/bash

#SBATCH -J bump_001         # job name
#SBATCH -o bump_001.o%j     # output and error file name (%j expands to jobID)
#SBATCH --ntasks=4
#SBATCH -N 2                 # number of nodes requested
#SBATCH --tasks-per-node=2   # Each task gets exactly one GPU
#SBATCH -t 00:10:00          # run time (hh:mm:ss) - 12.0 hours in this example.


# --------------------------------------------------------
# Might not be needed ?  
# --------------------------------------------------------
# #SBATCH --cpus-per-task=4    # Used for multithreading
# #SBATCH --gres=gpu:2         # Doesn't seem to work on Redhawk
# #SBATCH -n 2                 # total number of cpus requested.
# #SBATCH --exclusive          # request exclusive usage of your nodes. 
# ---------------------------------------------------------


# Execute the program:

# cd /home/donnacalhoun/projects/ForestClaw/code/forestclaw/applications/clawpack/shallow/2d/bump_cuda

# mpirun bump --user:cuda=T --clawpack46:mthlim="1" --clawpack46:order="2 2" --nout=100

#srun --mpi=pmix_v2 bump \
#     --user:cuda=T \
#     --cudaclaw:mthlim="4" \
#     --cudaclaw:order="2 2" \
#     --clawpack46:mthlim="0" \
#     --clawpack46:order="2 2"

nvprof -o bump_prof bump --user:cuda=T --cudaclaw:mthlim="4" --cudaclaw:order="2 2" 

