#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=3G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -t 1-24

#######################################################################

# Created by: Charley 2021, adapted from Chema
# Run an array to curate intron junctions with Portcullis, per bam file

#######################################################################

### Set up environment ###

module load anaconda3
conda activate portcullis

READS_DIR=/data/scratch/btx902/02_Gene_Evidences_RNA_Seq/Trimmed_Reads
ASS=/data/scratch/btx902/02_Gene_Evidences_RNA_Seq/Scaff/02_STAR_Alignments/CarCar_QM_v1_2021_12_Scaff_SoftMasked.fasta
BAM_DIR=/data/scratch/btx902/02_Gene_Evidences_RNA_Seq/Scaff/02_STAR_Alignments/bams

### Extract sample ID ###

Sample_File=$READS_DIR/ReadList.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)

### Run portcullis per bam file ###

portcullis full -t ${NSLOTS} -v \
--bam_filter --orientation FR --strandedness firststrand \
-o "portcullis_${SAMPLE}" \
$ASS $BAM_DIR/"${SAMPLE}_Aligned.sortedByCoord.out.bam"


