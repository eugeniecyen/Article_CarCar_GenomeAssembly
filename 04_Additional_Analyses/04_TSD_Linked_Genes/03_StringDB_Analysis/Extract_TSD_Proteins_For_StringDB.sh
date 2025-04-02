#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # 1 core sufficient
#$ -l h_rt=2:0:0  # 
#$ -l h_vmem=1G   # 
#$ -j y # combine std err and std out

#######################################################################

# Created by: James Gilbert (https://github.com/evoJDG), 2024

# Extract protein sequences of SD Genes from C mydas (as C mydas is in stringDB, but C caretta isn't yet. Based on protein functions so it's fine to use C mydas

# Need to extract the protein sequences to make the protein map
# Because the nucleotide sequences were extracted by coordinate (rather than sequence)
# the reading frames might not be perfect
# so need to use the sequence to find the corresponding protein ID and extract from protein fasta

# 1. Make blast db from C mydas cds reference
# 2. Blast SD gene seqs against db, add gene names
# 3. find the correct protein for each gene
# 4. extract the sequence using the protein ID

#######################################################################

############################
### Prepare Environment ####
############################

## load modules
module load anaconda3/2020.02
conda activate Orthofinder

## set directories
DIR=/data/SBCS-EizaguirreLab/James/SD_Flex/Submission-2
SYNDIR=${DIR}/01_Synteny
PROTDIR=${DIR}/03_ProtMap
METADIR=${DIR}/Genome_SD_Flex/00_Metadata
SCRIPTDIR=${DIR}/Genome_SD_Flex/01_Synteny/

# will save the green genome here:
CM_DIR=${SYNDIR}/data/Genomes/Cmydas

#######################
### Find Sequences ####
#######################
META=${METADIR}/Bentley_TSDgene_info.txt
GENES=$( cut -f 1 $META | tail -n +2 | sort | uniq )

# the reference assembly protein sequences
PROT_REF=${CM_DIR}/ncbi_dataset/data/GCF_015237465.2/protein.faa

CM_REF_CDS=${CM_DIR}/ncbi_dataset/data/GCF_015237465.2/cds_from_genomic.fna

# output file with protein IDs
PROTSEQIDS=${CM_DIR}/ncbi_dataset/data/GCF_015237465.2/SD_Prot_SeqIDs.txt
>$PROTSEQIDS

# look up each gene by name and print to output file
for GENE in $GENES; do

    # note the specific search and processing
    PROT_SEQ=$( grep -i "gene=$GENE]" $CM_REF_CDS \
    | head -n 1 | sed 's/.*protein_id=//g' \
    | cut -d ']' -f 1 )

    # if there are not matches i.e. prot_seq is empty, call it NA
    if [[ -z "$PROT_SEQ" ]]; then

        PROT_SEQ=NA

    fi

    echo -e "$GENE\t$PROT_SEQ" >> $PROTSEQIDS

done

##### add seqs for genes not found by name
## looking up seqs manually for these
# A2M
sed -i "/A2M/c\A2M\tXP_037768716.1" $PROTSEQIDS

# ANPEP
sed -i "/ANPEP/c\ANPEP\tXP_037766646.2" $PROTSEQIDS

# CA13
sed -i "/CA13/c\CA13\tXP_037748278.1" $PROTSEQIDS

# CIRBP
sed -i "/CIRBP/c\CIRBP\tXP_037741558.1" $PROTSEQIDS

# CNP
sed -i "/CNP/c\CNP\tXP_037742568.2" $PROTSEQIDS

# CPA3
sed -i "/CPA3/c\CPA3\tXP_007057453.2" $PROTSEQIDS

# CYP19A
sed -i "/CYP19A/c\CYP19A\tXP_037766849.1" $PROTSEQIDS

# DAX1
sed -i "/DAX1/c\DAX1\tXP_007063181.2" $PROTSEQIDS

# FAM132A
sed -i "/FAM132A/c\FAM132A\tXP_037737348.1" $PROTSEQIDS

# FAP
sed -i "/FAP/c\FAP\tXP_037769171.1" $PROTSEQIDS

# IFIT5
sed -i "/IFIT5/c\IFIT5\tXP_007071654.2" $PROTSEQIDS

# MPO
sed -i "/MPO/c\MPO\tXP_037736756.2" $PROTSEQIDS

# SCAP
sed -i "/SCAP/c\SCAP\tXP_037747260.1" $PROTSEQIDS

# SLC6A18
sed -i "/SLC6A18/c\SLC6A18\tXP_007052716.1" $PROTSEQIDS

# UCP2
sed -i "/UCP2/c\UCP2\tXP_037735488.1" $PROTSEQIDS


########## This produces the protein sequence ID for each gene
#### Now extract and print to file
PROTSEQS=${CM_DIR}/ncbi_dataset/data/GCF_015237465.2/SD_Prot_Sequences.faa
> $PROTSEQS

while read line; do

    # find the gene name and protein sequence ID
    GENE=$( echo $line | cut -d ' ' -f 1 )
    PROT_SEQ=$( echo $line | cut -d ' ' -f 2 )

    # first get only the protein sequence on one line
    SEQ=$( seqkit grep -p $PROT_SEQ $PROT_REF |  tail -n +2 | tr -d '\n' )

    # echo >gene name and seq to file
    echo -e ">$GENE\n${SEQ}" >> $PROTSEQS

done < $PROTSEQIDS














