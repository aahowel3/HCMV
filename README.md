# HCMV SLiM
on SOL, module load mamba, source activate HCMV 
have to install slim like every time even in HCMV environment 
/home/aahowel3/HCMV_slim 

/home/aahowel3/HCMV_slim/patient_data/patient_data.sh - general look of commands but mostly running line by line in align_quick.sh because its so buggy and breaks every line - but general idea in case we have to repeat it for other patients besdies B103

# Comparitive Pipeline
on Agave, plot testing playground in: /scratch/jemurra3/HCMV/trimmomatic/default_param/BWA/default 
source activate myenv

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
