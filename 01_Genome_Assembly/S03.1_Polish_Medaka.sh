#!/bin/bash
#$ -pe smp 12
#$ -l h_vmem=20G
#$ -l h_rt=240:00:00
#$ -l highmem
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Create polished consensus sequence with ONT's Medaka tool

#######################################################################

# Prep environment
module load anaconda3
conda activate medaka

READS_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/02_Quality_Control/03_Trimmed_Reads/02_NanoFilt
ASSEMBLY_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/04_Polished_Assembly/Racon/Wtdbg2/round_4
OUT_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/04_Polished_Assembly/Racon_Medaka/Wtdbg2/round_4

echo "Started: "`date`

# Run Medaka polishing
medaka_consensus \
-i $READS_DIR/turtle_porechop_q8minlen500.fastq.gz \
-d $ASSEMBLY_DIR/turtle_wtdbg2_racon_x4.fasta \
-o $OUT_DIR \
-m r941_prom_high_g4011 \
-t ${NSLOTS}

echo "Finished: "`date`
