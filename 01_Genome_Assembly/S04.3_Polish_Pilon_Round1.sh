#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=360G
#$ -l h_rt=240:00:00
#$ -l highmem
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Perform Round 1 of polishing with Illumina reads via Pilon 

#######################################################################

# Prep environment
module load java/11.0.2

PILON=/data/home/btx902/privatemodules/pilon-1.24.jar
ASSEMBLY_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/04_Polished_Assembly/Medaka
BAM_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/04_Polished_Assembly/Pilon/sorted_bams/dedup_bams
OUT_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/04_Polished_Assembly/Pilon/flye_medaka_pilon

echo "Started: "`date`

# Run Pilon
java -Xmx300G -jar $PILON \
--genome $ASSEMBLY_DIR/turtle_flye_medaka.fasta \
--bam $BAM_DIR/SL_063_flye_medaka.sort.dedup.bam \ 
--output turtle_flye_medaka_pilon \
--outdir OUT_DIR \
--changes
