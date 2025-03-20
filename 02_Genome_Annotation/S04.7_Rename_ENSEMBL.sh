#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=6G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Rename genes in ENSEMBL format with prefix "CarCarScaff"

#######################################################################

module load anaconda3
conda activate agat

agat_sp_manage_IDs.pl \
--gff CarCar_PASA_ThirdStep_AGAT_noSTOP.gff3 \
--ensembl \
--prefix CarCarScaff \
--type_dependent \
--tair \
-o CarCar_QM_v1.21.12_Sc.gff3

