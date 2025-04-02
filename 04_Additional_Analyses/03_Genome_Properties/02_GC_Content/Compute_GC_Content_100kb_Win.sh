#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=50G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y

#######################################################################

# Created by: Charley, Jun 2023, following workflow in Bentley et al. (2023) PNAS
# Calculate GC content for reference individual (SLK063)

#######################################################################

ASSEMBLY=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc_Chr0.fasta
OUT_DIR=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/04_GC_Content

### Create bed file from genome FASTA ###

module load python
source /data/SBCS-EizaguirreLab/Charley/environments/pyfaidx/bin/activate
PYFAIDX=/data/SBCS-EizaguirreLab/Charley/environments/pyfaidx/bin

python $PYFAIDX/faidx --transform bed \
$ASSEMBLY > CarCar_QM_v1.21.12_Sc_Chr0.bed

### Calculate GC content in sliding windows ###

# 100kb windows sliding by 20kb to match heterozygosity calculation 

# Make windows
deactivate
module purge
module load bedtools

bedtools makewindows \
-b CarCar_QM_v1.21.12_Sc_Chr0.bed \
-w 100000 -s 20000 > CarCar_QM_v1.21.12_Sc_Chr0_100kbwin_20kbstep.bed

# Calculate GC content per window
bedtools nuc \
-fi $ASSEMBLY \
-bed CarCar_QM_v1.21.12_Sc_Chr0_100kbwin_20kbstep.bed | cut -f1-5 \
> CarCar_GC_Content_100kbWin_20kbStep.txt

