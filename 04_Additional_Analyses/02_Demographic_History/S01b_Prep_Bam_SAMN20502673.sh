#!/bin/bash
#$ -pe smp 6
#$ -l h_vmem=1G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, Jul 2023
# Prep analysis-ready BAM files from raw Illumina reads for PSMC
# For Brazilian individual SAMN20502673 (Vilaca et al., 2021)

#######################################################################

OUT_DIR=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/03_PSMC/bams
ASSEMBLY=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc_Chr0.fasta

cd $OUT_DIR

### Download publicly available reads ###

module load sratools/2.10.8

prefetch SRR15328383
fasterq-dump --threads ${NSLOTS} --progress SRR15328383 # 74.3Gb, Illumina paired end reads 

### Quality trim reads for >Q30 ###

module purge
module load trimgalore

# Run TrimGalore
trim_galore --paired -q 30 --gzip --cores ${NSLOTS} \
SRR15328383_1.fastq SRR15328383_2.fastq

### Create BAM alignment files ###

module purge
module load bwa
module load samtools

# Align using BWA-MEM
# Convert SAM to sorted BAM
bwa mem -M -t ${NSLOTS} \
$ASSEMBLY \
$SRR15328383_1_val_1.fq.gz SRR15328383_2_val_2.fq.gz | \
samtools sort - > SRR15328383.sort.bam

# Index bam
samtools index SRR15328383.sort.bam > SRR15328383.sort.bam.bai

### Mark duplicates, calculate mean depth and mapping stats ###

module purge
module load java/11.0.2
PICARD=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/picard/build/libs/picard.jar

# Mark duplicates
java -jar $PICARD MarkDuplicates \
I=SRR15328383.sort.bam \
O=SRR15328383.sort.dedup.bam \
M=SRR15328383.marked_dup_metrics.txt \
ASO=coordinate

# Generate mean depth stats
java -jar $PICARD CollectWgsMetrics \
I=SRR15328383.sort.dedup.bam \
O=SRR15328383_wgsmetrics.txt \
R=$ASSEMBLY

# Re-index bam   
module purge
module load samtools

for i in *.dedup.bam
do
  samtools index $i > "${i}.bai"
done

# Generate mapping stats
for i in *.dedup.bam
do
  samtools flagstat $i > "${i}.flagstat"
done


