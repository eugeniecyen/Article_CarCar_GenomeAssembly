#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=3G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2024
# Run minimap2 with DGENIES suggested parameters to produce whole genome 
# alignments between our 'CarCar_QM_v1.21.12_Sc_Chr0' assembly and publicly 
# available chromosome-level assemblies for other sea turtle species: green
# and leatherback from Bentley et al. (2023); hawksbill from Guo et al., (2023)

#######################################################################

# Prep environment
module load anaconda3
conda activate minimap2

CARCAR_GEN=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc_Chr0.fasta
CHEMYD_GEN=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/01_Synteny_DGENIES/rCheMyd1.pri.v2_genomic.fna
DERCOR_GEN=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/01_Synteny_DGENIES/rDerCor1.pri.v4_genomic.fna
EREIMB_GEN=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/01_Synteny_DGENIES/GCA_030012505.1_ASM3001250v1_genomic.fna
OUT_DIR=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/01_Synteny_DGENIES/minimap

cd $OUT_DIR

# Loggerhead vs green
minimap2 -f 0.02 -t ${NSLOTS} \
$CARCAR_GEN $CHEMYD_GEN > $OUT_DIR/CarCar_CheMyd_dgenies_f0.02.paf

# Loggerhead vs leatherback
minimap2 -f 0.02 -t ${NSLOTS} \
$CARCAR_GEN $DERCOR_GEN > $OUT_DIR/CarCar_DerCor_dgenies_f0.02.paf

# Loggerhead vs hawksbill
minimap2 -f 0.02 -t ${NSLOTS} \
$CARCAR_GEN $EREIMB_GEN > $OUT_DIR/CarCar_EreImb_dgenies_f0.02.paf
