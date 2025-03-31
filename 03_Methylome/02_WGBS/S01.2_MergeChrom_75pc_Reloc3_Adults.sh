#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=26G
#$ -l h_rt=1:00:00
#$ -m bea

#######################################################################

# Created by Charley, 2023
# Bash submission script for R01.2_MergeChrom_75pc_Reloc3_Adults.R

#######################################################################

# Prep environment
module load R/4.2.2


# Call R script
Rscript ../R01.2_MergeChrom_75pc_Reloc3_Adults.R ${NSLOTS}
