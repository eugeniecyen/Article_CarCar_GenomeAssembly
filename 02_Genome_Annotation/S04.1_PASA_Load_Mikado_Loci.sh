#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=6G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Create PASA database for Mikado loci

#######################################################################

module load anaconda3
conda activate pasa

PATH=$PATH:/data/SBCS-EizaguirreLab/Charley/modules/conda/envs/pasa/opt/pasa-2.5.2

Launch_PASA_pipeline.pl \
-c pasa.alignAssembly.cfg \
-C -R -g CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta \
-t Clean_Mikado_Loci/mikado.transcripts.no.ncRNA.fa.clean \
-T -u Clean_Mikado_Loci/mikado.transcripts.no.ncRNA.fa \
--ALIGNERS blat --CPU ${NSLOTS}
