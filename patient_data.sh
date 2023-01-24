#!/bin/bash
#conda install -c bioconda entrez-direct
module load bwa/0.7.17
module load samtools/1.15.1
module load emboss/6.6.0


ACC=NC_006273.2
SRR=SRR1049475

# Shortcut to read names.
R1=${SRR}_1.fastq
R2=${SRR}_2.fastq
mkdir -p refs
REF=refs/$ACC.fa
BAM=$SRR.bam

# Download data
efetch -db=nuccore -format=fasta -id=$ACC | seqret -filter -sid $ACC > $REF
bwa index $REF

samtools faidx $REF
TAG="@RG\tID:$SRR\tSM:$SRR\tLB:$SRR"

####2 Align and generate a BAM file.
bwa mem -t 10 -R $TAG $REF $R1 $R2 | samtools view -b - > SRR.raw.bam
sambamba sort SRR.raw.bam
sambamba markdup SRR.raw.sorted.bam $BAM
samtools index $BAM
samtools stats $BAM | grep ^SN | cut -f 2- > stats/SRR.rmdup.sorted.bam.stats

####callinf variants
freebayes -f refs/NC_006273.2.fa --ploidy 1 --min-coverage 15 --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous --report-monomorphic -0 -i -X -u --region NC_00627$bcftools stats 475_sam_dupl_rm_no5.vcf > 475_sam_dupl_rm_no5.stats

vcffilter -f "QUAL > 0 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0" 475_sam_dupl_rm_no5.vcf > qual0_475_sam_dupl_rm_no5.vcf
bcftools stats qual0_475_sam_dupl_rm_no5.vcf > qual0_475_sam_dupl_rm_no5.stats
