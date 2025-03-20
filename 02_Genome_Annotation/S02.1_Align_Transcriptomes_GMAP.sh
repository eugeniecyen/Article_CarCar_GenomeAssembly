#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=2G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley 2021, adapted from Chema
# Align transcriptomes to genome assembly with GMAP

#######################################################################

### Prep environment ###
gmap_build=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/gmap-2021-12-17/bin/gmap_build
gmap=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/gmap-2021-12-17/bin/gmap

GMAP_DIR=/data/scratch/btx902/02_Gene_Evidences_RNA_Seq/Scaff/01_GMAP_Alignments
RNA_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/10_Annotation/02_Gene_Evidences_RNA_Seq/Transcriptomes

#######################################################################

### Build GMAP index ###

echo '______________Build GMAP index______________'

$gmap_build -D $GMAP_DIR \
-d CarCar_QM_v1_2021_12_Scaff_SoftMasked \
CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta


### Align transcripts to genome ###

echo '______________Align transcripts to HernFern2021______________'

$gmap -D $GMAP_DIR \
-d CarCar_QM_v1_2021_12_Scaff_SoftMasked \
-f 3 -n 0 -x 50 -B 4 -t ${NSLOTS} \
--gff3-add-separators=0 \
$RNA_DIR/GIBB01.fa > HernFern_2021_Vs_Genome_gmap.gff3


echo '______________Align transcripts to AusCha______________'

$gmap -D $GMAP_DIR \
-d CarCar_QM_v1_2021_12_Scaff_SoftMasked \
-f 3 -n 0 -x 50 -B 4 -t ${NSLOTS} \
--gff3-add-separators=0 \
$RNA_DIR/CarCar_AusCha_Trinity.fasta > Trinity_AusCha_Vs_Genome_gmap.gff3
