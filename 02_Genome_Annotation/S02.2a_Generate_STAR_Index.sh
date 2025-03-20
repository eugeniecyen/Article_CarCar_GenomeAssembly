#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=4G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Generate STAR index for soft-masked genome assembly

#######################################################################

STAR=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/STAR

$STAR --runMode genomeGenerate \
--runThreadN ${NSLOTS} \
--genomeDir /data/scratch/btx902/02_Gene_Evidences_RNA_Seq/Scaff/02_STAR_Alignments/CarCar_QM_v1_2021_12_Scaff_genome_STAR \
--genomeFastaFiles CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta

