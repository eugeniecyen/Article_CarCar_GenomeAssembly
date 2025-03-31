#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=2G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2023
# Reference-guided scaffolding of our contig-level assembly 'CarCar_QM_v1.21.12' 
# with RagTag against 28 chromosomes of the 'GSC_CCare_1.0' loggerhead assembly 
# released by Chang et al., (2023)

#######################################################################

# Extract chromosomes only from 'GSC_CCare_1.0' assembly
module load samtools

samtools faidx GCF_023653815.1_GSC_CCare_1.0_genomic.fna \
-r chroms.txt \
-o GCF_023653815.1_GSC_CCare_1.0_chrom_only.fna

module unload samtools

# Prep environment for RagTag
module load anaconda3
conda activate ragtag


# Run scaffolding with unplaced contigs
ragtag.py scaffold \
GCF_023653815.1_GSC_CCare_1.0_chrom_only.fna \
 CarCar_QM_v1.21.12.fasta \
-t ${NSLOTS}

# Run scaffolding with unplaced contigs concatenated into 'chr0' 
# with 100bp gap padding
ragtag.py scaffold \
GCF_023653815.1_GSC_CCare_1.0_chrom_only.fna \
CarCar_QM_v1.21.12.fasta \
-t ${NSLOTS} -C
