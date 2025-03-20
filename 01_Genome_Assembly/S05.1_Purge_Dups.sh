#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=2G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by Charley, 2021, adapted from Chema
# Haploidise genome assembly using Purge_Dups with Illumina reads by removing 
# haplotigs and contig overlaps based on read depth and sequence similarity 

#######################################################################

### Prep environment ### 
module load python/3.8.5

PURGE_DUPS=/data/home/btx902/privatemodules/purge_dups/bin
ASSEMBLY=/data/SBCS-EizaguirreLab/Turtle_Genome/04_Polished_Assembly/Pilon/flye_medaka_pilonx2/turtle_flye_medaka_pilonx2.fasta
ASSEMBLY_SPLIT=turtle_flye_medaka_pilonx2.split.fasta
SELF_ALN_GENOME=turtle_flye_medaka_pilonx2.split.self.paf.gz
BAM=/data/SBCS-EizaguirreLab/Turtle_Genome/06_BlobTools/bams/SL_063_flye_medaka_pilonx2.sort.dedup.bam

### Calculate read depth with ngscstat ###
echo 'NGSCSTAT _______________________________________________________________________'

$PURGE_DUPS/ngscstat $BAM

### Calculate cutoffs ###
echo 'CALCUTS _______________________________________________________________________'

$PURGE_DUPS/calcuts TX.stat > cutoffs 2> calcults.log

### Split FASTA file ###
echo 'SPLIT_FA ______________________________________________________________________'

$PURGE_DUPS/split_fa $ASSEMBLY > $ASSEMBLY_SPLIT

module unload python/3.8.5

### Perform minimap2 self alignment-alignment ###
echo 'MINIMAP2 ______________________________________________________________________'

module load anaconda3
conda activate minimap2

minimap2 -x asm5 -DP -t ${NSLOTS} \
$ASSEMBLY_SPLIT $ASSEMBLY_SPLIT | gzip -c - > $SELF_ALN_GENOME
 
conda deactivate
module unload anaconda3

### Purge dups ### 
echo 'PURGE_DUPS ____________________________________________________________________'

module load python/3.8.5

$PURGE_DUPS/purge_dups -2 -T cutoffs -c TX.base.cov \
$SELF_ALN_GENOME > dups.bed 2> purge_dups.log

### Get sequences ###
echo 'GET_SEQS ______________________________________________________________________'

$PURGE_DUPS/get_seqs dups.bed $ASSEMBLY > purged.fa 2> hap.fa
