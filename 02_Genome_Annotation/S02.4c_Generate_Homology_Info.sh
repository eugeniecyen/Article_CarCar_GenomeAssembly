#!/bin/bash
#$ -pe smp 3
#$ -l h_vmem=10G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley 2021, adapted from Chema

# Generate homology info for Mikado predicted transcripts via Diamond BLAST 
# against SwissProt database. Recommended before generating consensus gene/
# protein set, so Mikado can use it to improve the annotation.

#######################################################################

DIAMOND=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/diamond

FASTA=/data/scratch/btx902/02_Gene_Evidences_RNA_Seq/Scaff/04_Mikado/mikado_prepared.fasta

$DIAMOND blastx \
--query $FASTA \
--max-target-seqs 5 --sensitive --index-chunks 1 \
--threads ${NSLOTS} \
--db UniProt_SwissProt_2022_03_02.dmnd \
--evalue 1e-6 --outfmt 5 \
--out mikado_diamond_maxintron1.24Mb.xml

