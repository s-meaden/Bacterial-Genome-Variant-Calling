#!/bin/sh

# Mapping pipeline for mutation rate analysis. 61 PA14 genomes under different
# evolutionary treatments. Map to reference, then call variants to reference.
# Use ancestral genotype to screen SNPs by treatment. Reads already
# trimmed by MicrobesNG service.

# Using samtools mpileup instead of GATK HaplotypeCaller


# Map to reference w/ BWA-MEM

ref=Pseudomonas_aeruginosa_UCBPP-PA14_109.fna
bwa index $ref

# Folder for alignments:
mkdir 20170801_Meaden1/alignments/

# Loop through each sample (extract basename)
# (not using reads that lost a mate during trimming)
for sample in `ls 20170801_Meaden1/reads/*_1_trimmed.fastq.gz`
do
dir="20170801_Meaden1/alignments/"
base=$(basename $sample "_1_trimmed.fastq.gz")
bwa mem -t 16 $ref 20170801_Meaden1/reads/${base}_1_trimmed.fastq.gz 20170801_Meaden1/reads/${base}_2_trimmed.fastq.gz > ${dir}/${base}.sam
echo $base
done

# Convert SAM to BAM
for sample in 20170801_Meaden1/alignments/*.sam
do
echo $sample
describer=$(echo ${sample} | sed 's/.sam//')
echo $describer
samtools view -bT $ref $sample > ${describer}.uns.bam

# Sort BAM file
samtools sort ${describer}.uns.bam ${describer}

# Index bam file
samtools index ${describer}.bam

# Remove intermediates:
rm ${describer}.uns.bam

# Remove duplicate reads. Samtools and picard interchangeable for these steps.
# Using picard for compatability w/ GATK HaplotypeCaller (although both should work)

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/picard.jar MarkDuplicates INPUT=${describer}.bam OUTPUT=${describer}_dedup.bam METRICS_FILE=${describer}_metrics.txt
done

# Call SNPs w/ mpileup

samtools mpileup -uf $ref 20170801_Meaden1/alignments/*dedup.bam | bcftools view -bvcg > var.raw.bcf
bcftools view var.raw.bcf | /usr/share/samtools/vcfutils.pl varFilter -D100 > var.flt.vcf
