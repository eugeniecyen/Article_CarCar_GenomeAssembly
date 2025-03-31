#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -l centos 		# Run before Apocrita Rocky upgrade in Feb 2025

#######################################################################

# Created by: Alice Balard, 2023 (https://github.com/alicebalard)
# Convert methylation call bam to bed file and destrand across CpG sites

#######################################################################

# Prep environment
module load samtools

REF='/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc.fasta'
BAM=/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/bambeds/turtleONTmerged_2025.sorted.bam

#REF='/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc_Chr0.fasta'
#BAM=/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/bambeds/turtleONTmerged.sorted.bam

cd /data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/bambeds

# Run ModKit to generate destranded bed file per CpG for reference genome
modkit pileup $BAM \
turtleONTmerged_2025.sorted.bam_combineStrands.bed \
--cpg --combine-strands --ref $REF

# Prep bed file for loading into R methylKit
sed -e 's/ /\t/g' turtleONTmerged.sorted.bam_combineStrands.bed > turtleONTmerged.sorted.bam_combineStrands_tabdelim.bed