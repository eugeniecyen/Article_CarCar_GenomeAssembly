#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Update PASA database with Augustus transcripts

#######################################################################

module load anaconda3
conda activate pasa

PATH=$PATH:/data/SBCS-EizaguirreLab/Charley/modules/conda/envs/pasa/opt/pasa-2.5.2/scripts

Load_Current_Gene_Annotations.dbi \
-c pasa.alignAssembly.cfg \
-g CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta \
-P CarCar_BRAKER1_Combined_Renamed_EVM.gff3

