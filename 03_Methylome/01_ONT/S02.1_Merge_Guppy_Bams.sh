#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -l centos 		# Run before Apocrita Rocky upgrade in Feb 2025

#######################################################################

# Created by: Alice Balard, 2023 (https://github.com/alicebalard)
# Merge methylation call bam files outputted by Guppy

#######################################################################

# Prep environment
module load samtools

#file1=/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/Guppy_outdir/pass/bam_runid_a5072ac4a77030d5f831e8b84bf2f01cbc0acf46_0_0.bam
file1=/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/Guppy_outdir_2025/pass/bam_runid_a5072ac4a77030d5f831e8b84bf2f01cbc0acf46_0_0.bam

# Fix header issue
samtools view -H $file1 > /data/scratch/btx915/allSeq.sam

# Append sequences from all files into a single file
find /data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/Guppy_outdir_2025/pass/ \
 -name \*.bam -exec samtools view {} \; >> /data/scratch/btx915/allSeq.sam

# Convert to bam file
samtools view -bh /data/scratch/btx915/allSeq.sam > /data/scratch/btx915/turtleONTmerged_2025.bam

# Sort bam file
samtools sort /data/scratch/btx915/turtleONTmerged_2025.bam -o /data/scratch/btx915/turtleONTmerged_2025.sorted.bam

# Index bam file
samtools index /data/scratch/btx915/turtleONTmerged_2025.sorted.bam

# Get statistics
modkit summary /data/scratch/btx915/turtleONTmerged_2025.sorted.bam