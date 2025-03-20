#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=1G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Clean up annotation with AGAT toolkit

#######################################################################

module load anaconda3
conda activate agat

### Remove identical isoforms ###

agat_convert_sp_gxf2gxf.pl \
-g CarCar_PASA_ThirdStep.gff3 \
-o CarCar_PASA_ThirdStep_AGAT.gff3

### Remove any in-frame stop codons ###

gffread=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/gffread-0.12.7

$gffread \
-E CarCar_PASA_ThirdStep_AGAT.gff3 \
-g CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta \
-V -H -v \
-o CarCar_PASA_ThirdStep_AGAT_noSTOP.gff3

### Generate gene stats ###

agat_sq_stat_basic.pl \
-i CarCar_PASA_ThirdStep_AGAT_noSTOP.gff3 \
-g CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta

### Generate protein FASTA ###

$gffread \
-E CarCar_PASA_ThirdStep_AGAT_noSTOP.gff3 \
-g CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta \
-y CarCar_PASA_ThirdStep_AGAT_noSTOP.fasta

### Quick BUSCO ###

conda deactivate
conda activate busco

busco -c ${NSLOTS} -m prot \
-i CarCar_PASA_ThirdStep_AGAT_noSTOP.fasta \
-o busco_CarCar_PASA_ThirdStep_AGAT_noSTOP \
-l sauropsida_odb10

