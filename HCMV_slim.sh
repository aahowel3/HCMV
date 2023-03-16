#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -n 1            # number of "tasks" (default: 1 core per task)
#SBATCH -t 3-00:10:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition
#SBATCH --mem=16000mb
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)

#run on sols source activate HCMV 
#slim HCMV_burnin.slim

#python sc2-full_stats_forkAH.py burnin.100 23500 
#example python sc2-full_stats_forkAH.py 9000 23500
#Inputs are the file stem that comes before .output.fix and .output.ms (9000.output.fix and 9000.output.ms) and the genome length (23500)


#slim HCMV_expgrowth_psi.slim
#python sc2-full_stats_forkAH.py plasma.100 23500
#python sc2-full_stats_forkAH.py urine.100 23500


while IFS="," read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11
	do
	slim -d DFE="$col1" -d f0="$col2" -d f1="$col3" -d f2="$col4" -d f3="$col5" -d recombrate="$col6" -d recomb="$col7" -d murate="$col8" -d mu="$col9" -d progeny="$col10" -d gr="$col11" HCMV_burnin_final.slim
done < xaa

#should be params_burnin.txt but you're splitting it into xaa,xab etc with split --lines=72 params_burnin.txt

#while IFS="," read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12
#        do
#        slim -d DFE="$col1" -d f0="$col2" -d f1="$col3" -d f2="$col4" -d f3="$col5" -d recombrate="$col6" -d recomb="$col7" -d murate="$col8" -d mu="$col9" -d progeny="$col10" -d burnin="$col11" -d gr="$col12" HCMV_congenital_final.slim 
#done < params_congenital_tester.txt
