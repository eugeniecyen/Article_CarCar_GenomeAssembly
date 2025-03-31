#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -l centos 		# Run before Apocrita Rocky upgrade in Feb 2025

#######################################################################

# Created by: Alice Balard, 2023 (https://github.com/alicebalard)
# Count reads and output stats for BAM file

#######################################################################

# Prep environment
module load samtools

#BAM=/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/bambeds/turtleONTmerged.sorted.bam
BAM=/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/bambeds/turtleONTmerged_2025.sorted.bam

# Output BAM stats
echo "Number of reads:"

samtools view -c $BAM

echo "Total number of reads, number of mapped reads, number of properly paired reads (for paired-end data):"

samtools flagstat $BAM

echo "Mapping efficiency = (Number of mapped reads / Total number of reads) * 100"

echo "Calculate coverage depth:"

samtools depth $BAM | awk '{sum+=$3} END {print "Average coverage = ",sum/NR}'

samtools depth $BAM | cut -f3 | sort -n > temp2

median=$(( $(sed -n "$(($(wc -l < temp2)/2))p;$(($(wc -l < temp2)/2+1))p" temp2 | paste -sd+ | bc) / 2 ))

echo "The median is: $median"