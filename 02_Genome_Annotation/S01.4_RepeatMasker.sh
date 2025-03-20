#!/bin/bash
#$ -pe smp 20
#$ -l h_vmem=3G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley 2021, adapted from Chema
# Run RepeatMasker with custom repeat library to create soft masked genome

#######################################################################

### Prep environment ###

dfam_tetools=/data/home/btx902/singularity_containers/dfam-tetools.sh

ASS=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1_2021_12_Scaff.fasta
LIB=/data/scratch/btx902/01_Repeat_Masking/03_RepeatMasker/TEclass_renaming/CarCar-families-filt-TEclass.fa

# By default, RepeatMasker will start 2 threads for every slot requested,
# resulting in badly overloaded jobs. To avoid this, set a variable half of
# the allocated slots
REPCORES=$((NSLOTS / 2))

#######################################################################

### Run RepeatMasker ###

$dfam_tetools --singularity -- RepeatMasker \
-e rmblast -pa ${REPCORES} \
-xsmall -a -gff \
-lib $LIB $ASS

