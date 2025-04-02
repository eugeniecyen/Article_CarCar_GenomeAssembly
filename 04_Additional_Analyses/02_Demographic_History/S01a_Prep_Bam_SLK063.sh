#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=1G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, Jul 2023
# Prep analysis-ready BAM files from raw Illumina reads for PSMC
# For reference individual (SLK063)

#######################################################################

DATA_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/05_Illumina_Reads
OUT_DIR=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/03_PSMC/bams
ASSEMBLY=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc_Chr0.fasta

### Quality trim reads for >Q30 ###

cd $DATA_DIR

module load trimgalore

# Run TrimGalore
for i in *_1.fq.gz
do
  trim_galore --paired -q 30 --cores ${NSLOTS} \
  "${i:0:36}_1.fq.gz" "${i:0:36}_2.fq.gz"
done

cd $OUT_DIR

### Create BAM alignment files ###

module purge
module load bwa
module load samtools

# Align with BWA-MEM
# Convert SAM to sorted BAM
for i in $DATA_DIR/*_1.fq.gz

do
   bwa mem -M -t ${NSLOTS} \
   -R "@RG\tID:${i:81:9}.${i:92:1}\tSM:${i:57:6}\tPL:illumina\tPU:${i:81:9}.${i:92:1}.${i:57:6}" \
   $ASSEMBLY_DIR.fasta \
   $DATA_DIR/"${i:57:37}1_val_1.fq.gz" $DATA_DIR/"${i:57:37}2_val_2.fq.gz" | \
   samtools sort - > $OUT_DIR/"${i:57:6}_${i:81:9}.sort.bam"
done

# Index bams
for i in *.bam
do
   samtools index $i > "${i}.bai"
done

### Mark duplicates, merge BAMs, calculate mean depth and mapping stats ###

module purge
module load java/11.0.2
PICARD=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/picard/build/libs/picard.jar

# Mark duplicates and merge bams
java -jar $PICARD MarkDuplicates \
I=SL_063_H7L77DSX2.sort.bam \
I=SL_063_H7LJYDSX2.sort.bam \
O=SLK063.sort.dedup.bam \
M=SLK063.marked_dup_metrics.txt \
ASO=coordinate

# Generate mean depth stats
java -jar $PICARD CollectWgsMetrics \
I=SLK063.sort.dedup.bam \
O=SLK063_wgsmetrics.txt \
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





