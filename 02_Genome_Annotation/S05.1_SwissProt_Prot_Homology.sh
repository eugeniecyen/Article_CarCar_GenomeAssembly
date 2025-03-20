#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=2G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2022
# Assign functional annotation to genes via homology against curated SwissProt database

#######################################################################

DIR=/data/scratch/btx902/06_Functional_Annotation
DIAMOND=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/diamond

cd $DIR

$DIAMOND blastp \
-d UniProt_SwissProt_2022_03_02.dmnd \
-q CarCar_QM_v1_2021_12_Scaff_Annotation.aa.fa \
-o CarCar_QM_v1_2021_12_Scaff_Annotation_SwissProt.blastp \
--max-target-seqs 1 --evalue 1e-6 --max-hsps 1 -f 6 \
--threads ${NSLOTS}
