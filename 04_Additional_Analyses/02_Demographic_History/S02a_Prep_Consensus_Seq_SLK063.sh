#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=1G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -t 2-29	# Ignore chr0, as this consists of unplaced scaffold

#######################################################################

# Created by: Charley, Jul 2023
# Call consensus sequence in an array by chromosome
# For reference individual (SLK063)

#######################################################################

DATA_DIR=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/03_PSMC/bams
OUT_DIR=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/03_PSMC/01_Consensus_Sequences
ASSEMBLY=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc_Chr0.fasta

### Extract chrom IDs ###

Chrom_File=chrom_list.txt
CHROM=$(sed -n "${SGE_TASK_ID}p" $Chrom_File)

### Produce consensus sequence ### 

bcftools mpileup -Q 30 -q 30 -Ou \
-f $REF_GENOME -r $CHROM $DATA_DIR/SLK063.sort.dedup.bam | \
bcftools call -c | \
vcfutils.pl vcf2fq -Q 30 -d 17.2 -D 103.4 > $CHROM.fq

# Zip file
gzip $CHROM.fq
