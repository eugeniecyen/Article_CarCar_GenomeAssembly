#!/bin/bash
#$ -pe smp 8
#$ -l h_vmem=2G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley 2021, adapted from Chema

# Use Mikado serialise and pick to take all evidences (alignments, BLAST
# homologies, intron junctions, ORF predictions etc. ) to create consensus 
# set of best-fit gene models.

# NB. Fails when supplying full paths -> just use filenames 

#######################################################################

module load anaconda3
conda activate mikado

# Mikado serialise
mikado serialise \
--procs 1 \
--json-conf configuration_max_intron_1.24Mb.yaml \
--xml mikado_diamond_maxintron1.24Mb.xml \
--orfs mikado_prepared.fasta.transdecoder.bed \
--blast_targets UniProt_SwissProt_2022_03_02.fasta \
--transcripts mikado_prepared.fasta \
--junctions CarCar_junctions_consensus.bed

# Mikado pick
mikado pick \
--procs 10 \
--json-conf configuration_max_intron_1.24Mb.yaml
