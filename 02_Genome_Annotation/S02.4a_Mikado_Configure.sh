#!/bin/bash
#$ -pe smp 2
#$ -l h_vmem=1G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Generate Mikado configuration file

#######################################################################

module load anaconda3
conda activate mikado

mikado configure -t 2 \
--list mikado_config_input.tsv \
--reference CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta \
--mode permissive \
--scoring mammalian.yaml \
--copy-scoring mammalian.yaml \
--junctions CarCar_junctions_consensus.bed \
-bt UniProt_SwissProt_2022_03_02.fasta \
configuration.yaml
