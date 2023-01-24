
#!/bin/bash

#SBATCH -N 1  # number of nodes
#SBATCH -n 1  # number of "tasks" (default: allocates 1 core per task)
#SBATCH -t 7-00:00:00   # time in d-hh:mm:ss
#SBATCH --mem=16000mb
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --export=NONE   # Purge the job-submitting shell environment

module load picard/2.23.7

java -jar /packages/7x/picard/2.23.7/build/libs/picard.jar MarkDuplicates MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M=metrics.out I=PAV6_default_trimmomatic_BWA_default_consensus_sorted.bam O=PAV6_default_trimmomatic_BWA_default_consensus_sorted_rmmarkdups.bam
java -jar /packages/7x/picard/2.23.7/build/libs/picard.jar MarkDuplicates MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M=metrics2.out I=PAV6_default_trimmomatic_BWA_default_merlin_sorted.bam O=PAV6_default_trimmomatic_BWA_default_merlin_sorted_rmmarkdups.bam

