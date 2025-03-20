#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley 2021, adapted from Chema

# Use TransDecoder to predict ORFs and define the coding region for each
# transcript identified by Mikado prepare. This is recommended as many
# metrics used by Mikado to evaluate and rank transcripts rely on definition
# of the CDS.

#######################################################################

module load anaconda3
conda activate transdecoder

FASTA=/data/scratch/btx902/02_Gene_Evidences_RNA_Seq/Scaff/04_Mikado/mikado_prepared.fasta

# Step 1: Extract the long open reading frames
TransDecoder.LongOrfs -t $FASTA

# Step 2: Predict the likely coding regions
TransDecoder.Predict -t $FASTA

