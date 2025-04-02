#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=4G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Run RepeatModeler for genome assembly 'CarCar_QM_v1.21.'12_Sc'

#######################################################################

### Prep environment ###

dfam_tetools=/data/home/btx902/singularity_containers/dfam-tetools.sh
ASS=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1_2021_12_Scaff.fasta

#######################################################################

### Run RepeatModeler ###

# Format FASTA files for use with RepeatModeler
echo 'BuildDatabase ________________________'

$dfam_tetools --singularity -- BuildDatabase \
-name CarCar -engine ncbi $ASS

# Create species-specific database 
echo 'RepeatModeler _______________________'

$dfam_tetools --singularity -- RepeatModeler \
-LTRStruct -threads ${NSLOTS} -database CarCar

