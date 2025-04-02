#!/bin/bash
#$ -pe smp 6
#$ -l h_vmem=1G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, Nov 2022
# Extract Illumina read set enriched for mtDNA sequence by mapping against 
# published loggerhead mtDNA assembly (Drosopoulou et al., 2012), then
# assemble and annotate mitochondrial assembly with MitoZ

#######################################################################

DATA_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/05_Illumina_Reads
ENRICH_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/11_Mitogenome/01_Enrichment
MITO_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/11_Mitogenome/02_MitoZ_Run

### Extract reads enriched for mtDNA sequence ###

module load bwa
module load samtools

cd $ENRICH_DIR

# Align against mtDNA assembly using BWA-MEM
# Convert SAM to sorted BAM
bwa mem -M -t ${NSLOTS} \
CarCar_MT_Assembly_NC_016923.1.fa \
$DATA_DIR/SL_063_EDSW210006183-1a_H7LJYDSX2_L2_1.fq.gz \
$DATA_DIR/SL_063_EDSW210006183-1a_H7LJYDSX2_L2_2.fq.gz | \
samtools sort -n - > $ENRICH_DIRCarCar_MT_Assembly_NC_016923.1.qsort.bam

# Extract mapped reads only
# -F 4 → only include reads with none of the FLAGS in INT present
samtools view -b -F 4 -@ ${NSLOTS} \
CarCar_MT_Alignment_qsort.bam > CarCar_MT_Alignment_qsort_mapped.bam

# Convert BAM to FASTQ files
module unload bwa
module unload samtools
module load bedtools

bedtools bamtofastq \
-i CarCar_MT_Alignment_qsort_mapped.bam \
-fq CarCar_MT_Alignment_qsort_mapped_1.fq \
-fq2 CarCar_MT_Alignment_qsort_mapped_2.fq

gzip CarCar_MT_Alignment_qsort_mapped_1.fq
gzip CarCar_MT_Alignment_qsort_mapped_2.fq

### Assembly and annotate mt genome ###

cd $MITO_DIR

module unload bedtools
MITOZ=/data/home/btx902/singularity_containers/mitoz_3.4.sif # Set MitoZ Singularity package

singularity run $MITOZ mitoz all \
--thread_number ${NSLOTS} \
--outprefix CarCar_qsort_mito \
--clade Chordata --requiring_taxa Chordata \
--genetic_code 2 \
--fq1 $ENRICH_DIR/CarCar_MT_Alignment_qsort_mapped_1.fq.gz \
--fq2 $ENRICH_DIR/CarCar_MT_Alignment_qsort_mapped_2.fq.gz \
--assembler megahit \
--data_size_for_mt_assembly 0 # Set to 0 to use all data, NB. >8Gb not recommended









