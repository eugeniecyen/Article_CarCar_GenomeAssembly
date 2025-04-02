#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=50G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -t 2-29	# Ignore chr0, as this consists of unplaced scaffold

#######################################################################

# Created by: Charley, Jun 2023, following workflow in Bentley et al. (2023) PNAS
# Calculate genome-wide heterozygosity for reference individual (SLK063)
# Run GATK for variant calling with monomorphic sites in an array by chromosome

#######################################################################

BAM_DIR=/data/scratch/btx902/Heterozygosity/Bams # NB. Bams already produced during genome assembly 
GVCF_DIR=/data/scratch/btx902/Heterozygosity/Gvcfs
VCF_DIR=/data/scratch/btx902/Heterozygosity/Vcfs
ASSEMBLY=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc_Chr0.fasta

### Extract chrom IDs ###

Chrom_File=chrom_names.txt
CHROM=$(sed -n "${SGE_TASK_ID}p" $Chrom_File)

### Run GATK HaplotypeCaller ### 

gatk --java-options "-Xmx100G -Xms10G" HaplotypeCaller \
-I $BAM_DIR/SLK063.sort.dedup.bam \
-O $GVCF_DIR/"${CHROM}.mono.g.vcf.gz" \
-R $ASSEMBLY \
-ERC BP_RESOLUTION \
--output-mode EMIT_ALL_ACTIVE_SITES \
-L ${CHROM} -mbq 20

### Run GenotypeGVCFs with monomorphic sites ###

gatk --java-options "-Xms80G -Xmx100G" GenotypeGVCFs \
-R $ASSEMBLY \
-V $GVCF_DIR/${CHROM}.mono.g.vcf.gz \
-O $VCF_DIR/${CHROM}.mono.vcf.gz \
-L ${CHROM} \
--include-non-variant-sites \
--heterozygosity 0.0017 \
--standard-min-confidence-threshold-for-calling 0

### Remove unused alternate alleles ###

gatk --java-options "-Xms60G -Xmx90G" SelectVariants \
-R $ASSEMBLY \
-V $VCF_DIR/${CHROM}.mono.vcf.gz \
-O $VCF_DIR/${CHROM}.mono.TrimAlt.vcf.gz \
--remove-unused-alternates \
-L ${CHROM}





