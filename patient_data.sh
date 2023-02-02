#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -n 1            # number of "tasks" (default: 1 core per task)
#SBATCH -t 2-00:10:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition
#SBATCH --mem=10G
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)

#when running something on a high mem node use sbatch -p highmem scriptname.sh 
#and change --mem to 1002G 

ACC=NC_006273.2
SRR=SRR1049475

#Download fastq sequences 
#had to manually download SRR data - private ncbi servers you cant do fastq-dump 
#download it to a computer you dont care about its 4gigs
wget https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR1049475&display=download
gunzip SRR1049475.fastq.gz

# Download reference 
#will throw an error but still download it
efetch -db nuccore -format fasta -id $ACC > ${ACC}.fa
#nano NC_006273.2.fa FIRST to remove extra words in header or freebayes cant find chromo name
bwa index ${ACC}.fa
module load  samtools-1.16-gcc-11.2.0 #format required to load in SOL - ALWAYS STAYS UNMUTED
samtools faidx ${ACC}.fa



####2 Align and generate a BAM file. 
TAG="@RG\tID:$SRR\tSM:$SRR\tLB:$SRR" 
bwa mem -t 10 -R $TAG ${ACC}.fa -p ${SRR}.fastq > SRR.raw.sam
samtools view -S -b SRR.raw.sam > SRR.raw.bam
sambamba sort SRR.raw.bam
sambamba markdup SRR.raw.sorted.bam SRR.raw.sorted.mkdup.bam
samtools index SRR.raw.sorted.mkdup.bam
samtools stats SRR.raw.sorted.mkdup.bam | grep ^SN | cut -f 2- > SRR.raw.sorted.mkdup.bam.stats

####callinf variants
freebayes -f NC_006273.2.fa --ploidy 1 --min-coverage 15 --haplotype-length 0 --min-alternate-count 1 --min-alternate-fraction 0 --pooled-continuous --report-monomorphic -0 -i -X -u --region NC_006273.2:62500-86000 SRR.raw.sorted.mkdup.bam > 475_sam_dupl_rm_no5.vcf
bgzip 475_sam_dupl_rm_no5.vcf
tabix 475_sam_dupl_rm_no5.vcf.gz
bcftools stats 475_sam_dupl_rm_no5.vcf.gz > 475_sam_dupl_rm_no5.stats 
#mamba install -c bioconda vcflib
vcffilter -f "QUAL > 0 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0" 475_sam_dupl_rm_no5.vcf.gz > qual0_475_sam_dupl_rm_no5.vcf
bgzip qual0_475_sam_dupl_rm_no5.vcf
tabix qual0_475_sam_dupl_rm_no5.vcf.gz
bcftools stats qual0_475_sam_dupl_rm_no5.vcf.gz > qual0_475_sam_dupl_rm_no5.stats





####IGNORE JUST LEAVE IT INTERLEAVED AND USE BWA MEM -P
#headers are completely screwed - edit to remove extra characters so reformat.sh will work properly
#cat SRR1049475.fastq | sed 's/\s.*$//' > SRR1049475_renamed.fastq
#reformat.sh in=SRR1049475_renamed.fastq out1=SRR1049475_R1.fq out2=SRR1049475_R2.fq
#reformat.sh in=SRR1049475.fastq allowidenticalnames=t out1=SRR1049475_R1_renamed.fq out2=SRR1049475_R2_renamed.fq
#repair.sh -Xmx1000g overwrite=t in=SRR1049475_R1.fq in2=SRR1049475_R2.fq out=SRR1049475_R1_sorted.fq out2=SRR1049475_R2_sorted.fq outs=SRR1049475_singletons.fq
#############################################
