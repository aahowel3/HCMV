# HCMV SLiM
on SOL, source activate HCMV 

# Comparitive Pipeline
on Agave, plot testing playground in: /scratch/jemurra3/HCMV/trimmomatic/default_param/BWA/default 

markdups_tester.sh - create deduplicated .bams for additional coverage plots
also creates a DT tagged bam without removing them to categorize duplicate types

one liner I used to pull out the start position and CIGAR string from a bam and the Rcode to graph it by duplication number
samtools view PAV6_default_trimmomatic_BWA_default_consensus_sorted.bam | awk '
BEGIN{OFS=":"} {print $4, $6}' > PAV6_default_trimmomatic_BWA_default_consensus_sorted.dupstags.txt

pull out read name and DT tag
samtools view PAV6_default_trimmomatic_BWA_default_merlin_sorted_markdups.bam | awk '{for(i=1;i<=NF;i++){if($i~/^DT/){a=$i}} print $1,a}' | head
