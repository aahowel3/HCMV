# HCMV SLiM
on SOL, module load mamba, source activate HCMV 
PATIENT DATA 
/home/aahowel3/HCMV_slim/patient_data/patient_data.sh - general look of commands but mostly running line by line because its breaks - but general idea in case we have to repeat it for other patients besides B103

moved rerun to scratch bc I ran out of room on /home 
/scratch/aahowel3/HCMV_slim/patient_data_rerun
bwa -p does not work still reads as SE reads - to remove everything after first space and replace last . wiht a / use: head SRR1049475.fastq | sed 's/\s.*$//' | sed 's/\(.*\)\./\1\//'
AFTER I FIXED THE PAIRED END ERROR AND RAN IT THROUGH THE EXACT PIPELINE THE SUMMARY STAT SCRIPT FAILS NOW 
ALREADY TRIED REMOVING MULTIALLELIC SITES - issue with sites not all being the same coverage - fixed in R script HCMV_vcffiltering.R - have to downsample to 100 and remove 2% reads

new patient data HANCHILD 1,2,3,4 in /scratch/aahowel3/HCMV_slim/patient_data/hanchild

SIMUALTED DATA 
/home/aahowel3/HCMV_slim/simulations 
HCMV_slim.sh - runs burnin or expgrowth_psi through bash 
python sc2-full_stats_forkAH.py 9000 23500
Inputs are the file stem that comes before .output.fix and .output.ms (9000.output.fix and 9000.output.ms) and the genome length (23500)

on local - plot_summarystats.R for plotting 
on local - HCMV_vcffiltering.R - lets see if we can do the 100x downsampling and 2% removal reading in the vcf and spitting out a vcf - so that it goes cleanly back into terbots summary stats script

Patient data: SOL: /scratch/aahowel3/HCMV_slim/patient_data_rerun
Simulations: SOL: /home/aahowel3/HCMV_slim/simulations

Final simulations that iterate through all possible param combos in: [aahowel3@login02:~/HCMV_slim/simulations/final_sims]$ /home/aahowel3/HCMV_slim/simulations/final_sims

on local var_mu_recomb_rate.R generates the variable mutation and recombination rate files and generates the all possible param combinations params.txt file

# Comparitive Pipeline
on Agave, plot testing playground in: /scratch/jemurra3/HCMV/trimmomatic/default_param/BWA/default 
source activate myenv

on local - coverages_HCMVdata.R for plotting 

markdups_tester.sh - create deduplicated .bams for additional coverage plots
also creates a DT tagged bam without removing them to categorize duplicate types - creates a rmmarkdups.bam and markdups.bam

one liner I used to pull out the start position and CIGAR string from a bam and the Rcode to graph it by duplication number
samtools view PAV6_default_trimmomatic_BWA_default_consensus_sorted.bam | awk '
BEGIN{OFS=":"} {print $4, $6}' > PAV6_default_trimmomatic_BWA_default_consensus_sorted.dupstags.txt

pull out read name and DT tag
now in the script pull_duplicates since agave quits on you if you run on login node 
samtools view PAV6_default_trimmomatic_BWA_default_merlin_sorted_markdups.bam | awk '{for(i=1;i<=NF;i++){if($i~/^DT/){a=$i}} print $1,a}' | head

from pulled out read name and DT tag pull out those names from the fastqs to get shortened fastq that only came from dup reads
/scratch/jemurra3/HCMV/trimmomatic/default_param/BWA/default/crossref.sh

from those desired fastqs pull out first line each read to read into R can seperate by column on read-in
#example head ERR3013919_2_fastp.fastq | sed -n '1~4p'> ERR3013919_2_fastp.loc 
(myenv) [aahowel3@agave3:/scratch/jemurra3/HCMV/trimmomatic/default_param/BWA/default]$ sed -n '1~4p' out.fq > out_locs.txt

checking if any tag differences between R1 and R2 - getting more dup tags than possible tags in file
awk -F" " 'NR==FNR{a[FNR]=$1; next} {print $1, $1 == a[FNR] ? "ok" : "error"}' ERR_1.loc.txt ERR_2.loc.txt > check.out 
^^false alarm this was because you were looking at fastqs trimmed by a different tool than the pipeline that generated the alignment you were looking at 
