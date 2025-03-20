#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=1G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Round 1 of PASA annotation loading, comparison and updates to incorporate 
# transcript alignnments into gene structures

#######################################################################

module load anaconda3
conda activate pasa

PATH=$PATH:/data/SBCS-EizaguirreLab/Charley/modules/conda/envs/pasa/opt/pasa-2.5.2

Launch_PASA_pipeline.pl \
-c pasa.annotationCompare.cfg -A \
-g CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta \
-t Clean_Mikado_Loci/mikado.transcripts.no.ncRNA.fa.clean \
--CPU ${NSLOTS}
