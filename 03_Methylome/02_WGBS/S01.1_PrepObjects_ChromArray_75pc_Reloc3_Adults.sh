#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_vmem=2G
#$ -l h_rt=240:00:00
#$ -t 2-29


#######################################################################

# Created by Charley, 2023
# Bash submission script for R01.1_PrepObjects_ChromArray_75pc_Reloc3_Adults.R

#######################################################################

# Prep environment
module load R/4.2.2

# Call R script

# With 2 arguments to pass into R:
# 1. Chromosome number (which is the job array number i.e. SGE_TASK_ID). This must be the 1st arg.
# Note we have chromosomes 0:28 but 0 can't be an array ID so the array is 1:29, 
# then in the R script the chr number is set to SGE_TASK_ID - 1
# 2. Number of cores available for script ($NSLOTS)

Rscript ../R01.1_PrepObjects_ChromArray_75pc_Reloc3_Adults.R ${SGE_TASK_ID} ${NSLOTS}
