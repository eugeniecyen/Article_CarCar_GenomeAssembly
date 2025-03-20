#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 20
#$ -l h_vmem=1G
#$ -l h_rt=240:00:00
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Map trimmed Illumina reads to the genome assembly (Flye + Medaka)
# with BWA, index, mark duplicates and merge BAMs with Picard.

#######################################################################


### Prep BAM alignments from Illumina reads ###

# Prep environment
module load bwa
module load samtools
PICARD=/data/home/btx902/privatemodules/picard/build/libs/picard.jar

DATA_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/05_Illumina_Reads
ASSEMBLY_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/04_Polished_Assembly/Medaka
BAM_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/04_Polished_Assembly/Pilon/sorted_bams
DEDUP_DIR=$BAM_DIR/dedup_bams

# Align using BWA-MEM
# Convert SAM to sorted BAM
for i in $DATA_DIR/*.fq.gz

do
   bwa mem -M -t ${NSLOTS} \
   -R "@RG\tID:${i:81:9}.${i:92:1}\tSM:${i:57:6}\tPL:illumina\tPU:${i:81:9}.${i:92:1}.${i:57:6}" \
   $ASSEMBLY_DIR/turtle_flye_medaka.fasta \
   $DATA_DIR/"${i:57:37}1.fq.gz" $DATA_DIR/"${i:57:37}2.fq.gz" | samtools sort - > $BAM_DIR/"${i:57:6}_${i:81:9}_flye_medaka.sort.bam"
done

echo "Finished: "`date`

# Index BAMs
for i in *.bam
do
   samtools index $i > "${i}.bai"
done

### Mark duplicates and merge BAMs ###

# Prep environment
module unload bwa
module load java/11.0.2

# Run Picard MarkDuplicates 
java -jar $PICARD MarkDuplicates \
I=$BAM_DIR/SL_063_H7L77DSX2_flye_medaka.sort.bam \
I=$BAM_DIR/SL_063_H7LJYDSX2_flye_medaka.sort.bam \
O=$DEDUP_DIR/SL_063_flye_medaka_pilon.sort.dedup.bam \
M=$DEDUP_DIR/SL_063_flye_medaka_pilon.marked_dup_metrics.txt \
ASO=coordinate

# Generate mapping quality stats
for i in *.bam
do
   samtools flagstat $i > "${i}.flagstat"
done
