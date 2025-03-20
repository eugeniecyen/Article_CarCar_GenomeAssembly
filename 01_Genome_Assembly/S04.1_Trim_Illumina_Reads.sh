#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=6G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Remove adapters and low quality reads from Illumina data

#######################################################################

# Prep environment
module load trimgalore # v.0.6.5

DATA_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/05_Illumina_Reads

# Run TrimGalore
for i in $DATA_DIR/*.fq.gz
do
  trim_galore $i
done

echo "Finished: "`date`
