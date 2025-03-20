#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley 2021, adapted from Chema

# Generate a sorted, non-redundant GTF with all of the input assemblies, 
# using config file created in previous step, with max intron length 1.24 Mbp
# as identified in VGP CheMyd annotation

#######################################################################

module load anaconda3
conda activate mikado

mikado prepare --json-conf configuration_max_intron_1.24Mb.yaml
