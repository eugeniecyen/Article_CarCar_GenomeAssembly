#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=4G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Run BRAKER1 to perform gene prediction with RNA-Seq hints

#######################################################################

module load anaconda3
conda activate braker2

braker.pl \
--GENEMARK_PATH=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/braker_dependencies/ProtHint-2.6.0/dependencies/GeneMarkES \
--skipAllTraining --species=chicken \
--genome=CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta \
--bam=merged.bam \
--softmasking \
--verbosity=3 \
--cores ${NSLOTS}

# Run BUSCO 
cd braker

conda deactivate
conda activate busco

busco -c ${NSLOTS} -m prot \
-i augustus.hints.aa \
-o busco \
-l sauropsida_odb10

