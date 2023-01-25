# HCMV SLiM
on SOL, source activate HCMV 

# Comparitive Pipeline
on Agave, plot testing playground in: /scratch/jemurra3/HCMV/trimmomatic/default_param/BWA/default 

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
head ERR3013919_2_fastp.fastq | sed -n '1~4p'> ERR3013919_2_fastp.loc 
