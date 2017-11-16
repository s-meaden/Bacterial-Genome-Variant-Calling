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
#for sample in 20170801_Meaden1/alignments/*.sam
#do
#echo $sample
#describer=$(echo ${sample} | sed 's/.sam//')
#echo $describer
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
#done

# Change headers so compatible w/ HaplotypeCaller
#for i in 20170801_Meaden1/alignments/*dedup.bam
#do
#~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/picard.jar AddOrReplaceReadGroups I=${i} O=${i}_renamed.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${i}

# Index new BAM files:
#~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/picard.jar BuildBamIndex INPUT=${i}_renamed.bam
#done

# Create sequence dictionary:
#~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/picard.jar CreateSequenceDictionary R=$ref O=${ref}.dict

# Call SNPs w/ HaplotypeCaller. Make sure list of BAM files includes path
#~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar \
#  -T HaplotypeCaller \
#  -R $ref \
#  -I bam_list_hc.list \
#  -nct 16 \
#  --genotyping_mode DISCOVERY \
#  -newQual \
#  -stand_call_conf 10 \
#  -o raw_variants_hc.vcf \
#  --sample_ploidy 1

# Take raw variants and apply BSQR recalibraton using the SNPs w/ highest quality.
~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T SelectVariants -R $ref -V raw_variants_hc.vcf -selectType SNP -o raw_snps.vcf

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T SelectVariants -R $ref -V raw_variants_hc.vcf -selectType INDEL -o raw_indels.vcf

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o filtered_snps.vcf

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o filtered_indels.vcf

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T BaseRecalibrator -R $ref -I bam_list_hc.list -knownSites filtered_snps.vcf -knownSites filtered_indels.vcf -o recal_data.table

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T BaseRecalibrator -R $ref -I bam_list_hc.list -knownSites filtered_snps.vcf -knownSites filtered_indels.vcf -BQSR recal_data.table -o post_recal_data.table

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $ref -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T PrintReads -R $ref -I bam_list_hc.list -BQSR recal_data.table -o recal_reads.bam

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref -I bam_list_hc.list -o raw_variants_recal.vcf

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T SelectVariants -R $ref -V raw_variants_recal.vcf -selectType SNP -o raw_snps_recal.vcf

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T SelectVariants -R $ref -V raw_variants_recal.vcf -selectType INDEL -o raw_indels_recal.vcf

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V raw_snps_recal.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o filtered_snps_final.vcf

~/programs/jdk1.8.0_151/jre/bin/java -jar ~/programs/GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V raw_indels_recal.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o filtered_indels_recal.vcf

#java -jar snpEff.jar -v snpeff_db filtered_snps_final.vcf > filtered_snps_final.ann.vcf
