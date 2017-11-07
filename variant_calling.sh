#!/bin/sh

# Mapping pipeline for mutation rate analysis. 61 PA14 genomes under different
# evolutionary treatments. Map to reference, then call variants to reference.
# Use ancestral genotype to screen SNPs by treatment. Reads already
# trimmed by MicrobesNG service.


# Map to reference w/ BWA-MEM

ref=Pseudomonas_aeruginosa_UCBPP-PA14_109.fna
#bwa index $ref

# Folder for alignments:
#mkdir 20170801_Meaden1/alignments/

# Loop through each sample (extract basename)
# (not using reads that lost a mate during trimming)
#for sample in `ls 20170801_Meaden1/reads/*_1_trimmed.fastq.gz`
#do
#dir="20170801_Meaden1/alignments/"
#base=$(basename $sample "_1_trimmed.fastq.gz")
#bwa mem -t 16 $ref 20170801_Meaden1/reads/${base}_1_trimmed.fastq.gz 20170801_Meaden1/reads/${base}_2_trimmed.fastq.gz > ${dir}/${base}.sam
#echo $base
#done

# Convert SAM to BAM
for sample in 20170801_Meaden1/alignments/*.sam
do
echo $sample
describer=$(echo ${sample} | sed 's/.sam//')
echo $describer
#samtools view -bT $ref $sample > ${describer}.uns.bam

# Sort BAM file
#samtools sort ${describer}.uns.bam ${describer}

# Index bam file
#samtools index ${describer}.bam

# Remove intermediates:
#rm ${describer}.uns.bam

# Remove duplicate reads. Samtools and picard interchangeable for these steps.
# Using picard for compatability w/ GATK HaplotypeCaller (although both should work)

#~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/picard.jar MarkDuplicates INPUT=${describer}.bam OUTPUT=${describer}_dedup.bam METRICS_FILE=${describer}_metrics.txt
done

# Change headers so compatible w/ HaplotypeCaller
#for i in 20170801_Meaden1/alignments/*dedup.bam
#do
#~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/picard.jar AddOrReplaceReadGroups I=${i} O=${i}_renamed.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${i}

# Index new BAM files:
#~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/picard.jar BuildBamIndex INPUT=${i}_renamed.bam
#done

# Create sequence dictionary:
#~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/picard.jar CreateSequenceDictionary R=$ref O=${ref}.dict

# Call SNPs w/ HaplotypeCaller. Need to input all BAM files separately?
~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar \
  -T HaplotypeCaller \
  -R $ref \
  -I bam_list_hc.list \
  -nct 16 \
  --genotyping_mode DISCOVERY \
  -newQual \
  -stand_call_conf 10 \
  -o raw_variants_hc.vcf \
  --sample_ploidy 1

# Take raw variants and apply BSQR recalibraton using the SNPs w/ highest quality.
