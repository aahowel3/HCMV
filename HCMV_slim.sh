#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -n 1            # number of "tasks" (default: 1 core per task)
#SBATCH -t 2-00:10:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition
#SBATCH --mem=16000mb
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)

#module load slim/3.7.1
#slim HCMV_burnin.slim

#python sc2-full_stats_forkAH.py burnin.100 23500 
#example python sc2-full_stats_forkAH.py 9000 23500
#Inputs are the file stem that comes before .output.fix and .output.ms (9000.output.fix and 9000.output.ms) and the genome length (23500)

#slim HCMV_expgrowth_psi.slim
python sc2-full_stats_forkAH.py plasma.100 23500
python sc2-full_stats_forkAH.py urine.100 23500
