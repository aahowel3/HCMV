# HCMV SLiM

n each replicate how to get a "goodcombos" list to filter into terbotfiltering to fitler into sc2
cat 100bi_10tricombos.csv | grep "plasma_100_bi" | grep "good" | awk '{print $2}'

new file HCMV_slim_onlysim_leftovers2.sh - checks whichs combos are missing and runs only those
input file HAS to have carriage returns removed first
tr '\r' , < params_congenital.txt > params_congenital_nor.txt

# for combining all the final files to then put into nice plots
for x in simulations_replicates/replicate_*/good_combos/*.1000.filter.output.csv; do sed -n 1630p $x >> allfinalsnps_andcombos.csv; done
for x in simulations_replicates/replicate_*/good_combos/*.100.filter.output.csv; do sed -n 1626p $x >> allfinalsnps_andcombos.csv; done

where are our 1x replicates? lab computer #1 - in Documents/simulations/rescale ls *1.0.w*.fix.txt | wc -l - should be about 141 - need like 3 more 

on SOL, module load mamba, source activate HCMV 
PATIENT DATA 
/home/aahowel3/HCMV_slim/patient_data/patient_data.sh - general look of commands but mostly running line by line because its breaks - but general idea in case we have to repeat it for other patients besides B103

moved rerun to scratch bc I ran out of room on /home 
/scratch/aahowel3/HCMV_slim/patient_data_rerun
bwa -p does not work still reads as SE reads - to remove everything after first space and replace last . wiht a / use: head SRR1049475.fastq | sed 's/\s.*$//' | sed 's/\(.*\)\./\1\//'
AFTER I FIXED THE PAIRED END ERROR AND RAN IT THROUGH THE EXACT PIPELINE THE SUMMARY STAT SCRIPT FAILS NOW 
ALREADY TRIED REMOVING MULTIALLELIC SITES - issue with sites not all being the same coverage - fixed in R script HCMV_vcffiltering.R - have to downsample to 100 and remove 2% reads

new patient data HANCHILD 1,2,3,4 in /scratch/aahowel3/HCMV_slim/patient_data/hanchild

pipeline for B103 data is:
1. fastqs go through patientdata.sh and output is qual_nomulti.vcf
2. get that vcf into vcf HCMV_vcffiltering.R to downsample and 2% filter output is good.csv
3. That vcf goes into pythonsc2 and output is good.csv 
4. csv goes into plotsummarystats.R

Started running things locally since 1x gr sims take so long - pfeiferlab@10.210.91.237
Run sc2.py on ALL combos first
Ran a quick and dirty check on which param combos passed by doing: sed -n 1630p congenital_urine.DFE3.2.0e-07.9.8e-07.0.1.0.38.without.1000.output.csv
*****1630 for 1000x 1626 for 100x*****
Created a check_paramcombos.csv by the sed -n method (double check if plasma 1000, plasma 100, urine 1000, urine 100 lenghts are uneven - because the .sc2.py method will not work if there are ZERO segregating sites period 
Get a list of acceptable param combos using the chunk of code labelled "after filtering 2%" in HCMV_vcffiltering.R 

ONCE YOU GET THAT LIST
then in /Users/pfeiferlab/Documents/simulations/rescale_working/filter use copy_edit.sh to run acceptable combos through HCMV_filtering_forterbotscriptinput.R and sc2-full2_stats-v0.2_1_for100.py and plots_summarystats_updated_args.R for full results on filtered data. 

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

WHEN RUNNING FOR A VCF HAVE TO PUT ABSOLUTE LENGTH (87000) NOT RELATIVE LENGTH (23500)
python3 sc2-full_stats-v0.2.py qual0_475_sam_dupl_rm_no5 v 87000

THIS ONE WORKS
(base) pfeiferlab@tCLAsol-F5KTT0T1F694 simulations % python sc2-full2_stats-v0.2_1.py congenital_urine.DFE3.2.0e-06.0.0.01.0.38.without.1000 m 23500

get line that matters: sed -n 1630p  congenital_urine.DFE3.2.0e-07.9.8e-07.0.1.0.38.without.1000.output.csv
*****1630 for 1000x 1626 for 100x

squish all relavant pdfs together:  while read -r line; do gs -1 -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=${line}.plots.all.pdf congenital_plasma.${line}.100.plots.pdf congenital_plasma.${line}.1000.plots.pdf congenital_urine.${line}.100.plots.pdf congenital_urine.${line}.1000.plots.pdf; done < goodcombos2.txt

# what to do if runs dont finish and every replicate is missing a different combination
for each replicate get a file with replicate name and list of complete combos
for i in replic*; do cd $i; echo $i > ../checkfile_${i}.txt; ls congenital_urine*.100.output.fix >> ../checkfile_${i}.txt; cd ../; done
paste all those rep files into a columnwise rep file
paste -d' ' checkfile_replicate_{1..25}.txt > all_check.out
theres a section in mu_var.recomb.R script that you can cross check which are missing from which replicate

HCMV_vcffiltering2.R - this was your original HCMV_filtering script to get number of seg sites from empirical data and convert to vcf format terbots script would take (downsampled) - messed it up a bit with script at the top which has been moved to HCMV_filtering_refined - bottom part was most useful - dont know if youll need to convert downsampled emprical to anymore since youre not evaluating simualtions on summary statistics - only seg sites

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
