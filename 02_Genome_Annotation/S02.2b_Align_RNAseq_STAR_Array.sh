#!/bin/bash
#$ -j y
#$ -cwd
#$ -pe smp 10
#$ -l h_vmem=3G
#$ -l h_rt=24:00:00
#$ -t 1-24

#######################################################################

# Created by: Charley 2021, adapted from Chema
# Run an array aligning all RNAseq reads to genome with STAR 

#######################################################################

### Set up environment ###

STAR=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/STAR
stringtie=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/stringtie-2.2.1.Linux_x86_64/stringtie

GENOME_DIR=/data/scratch/btx902/02_Gene_Evidences_RNA_Seq/Scaff/02_STAR_Alignments/CarCar_QM_v1_2021_12_Scaff_genome_STAR
READS_DIR=/data/scratch/btx902/02_Gene_Evidences_RNA_Seq/Trimmed_Reads

### Extract sample ID ###

Sample_File=$READS_DIR/ReadList.txt
SAMPLE=$(sed -n "${SGE_TASK_ID}p" $Sample_File)

#######################################################################

### Align and convert to GTF ###
$STAR --runThreadN ${NSLOTS} \
--genomeDir $GENOME_DIR \
--readFilesCommand zcat \
--readFilesIn $READS_DIR/"${SAMPLE}_trim_1_paired.fastq.gz" $READS_DIR/"${SAMPLE}_trim_2_paired.fastq.gz" \
--outSAMtype BAM SortedByCoordinate \
--outWigType wiggle --outWigStrand Stranded \
--outFileNamePrefix "${SAMPLE}_"

$stringtie "${SAMPLE}_Aligned.sortedByCoord.out.bam" \
-o "${SAMPLE}_transcripts.gtf" \
-p ${NSLOTS} -v

