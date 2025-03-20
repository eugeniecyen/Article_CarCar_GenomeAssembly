#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley 2021, adapted from Chema
# Merge curated intron junction bed files with junctools

#######################################################################

### Set up environment ###

module load anaconda3
conda activate portcullis

PORTC_DIR=/data/scratch/btx902/02_Gene_Evidences_RNA_Seq/Scaff/03_Portcullis
cd $PORTC_DIR

bed_list=$(for i in portcullis_*/3-filt/portcullis*.bed; do echo "$i"; done)

### Run junctools set ###

junctools set -m 1 --operator max \
-o CarCar_junctions_consensus.bed \
consensus $bed_list
