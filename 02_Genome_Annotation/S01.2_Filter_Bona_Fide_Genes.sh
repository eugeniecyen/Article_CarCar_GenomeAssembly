#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=2G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2023
# Filter out bona fide genes from repeat library

# Using a curated VGP genome proteome for Chelonia mydas, as it is 
# the most complete annotation available for a sea turtle species 

#######################################################################

### Prep environment ###

diamond=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/diamond

#######################################################################

### BLAST for potential genes in our CarCar repeat library ###

$diamond blastx \
-d VGP_CheMyd_Prot_noTEs.dmnd \
-q CarCar-families.fa \
-o CarCar_VGPCheMyd_NoTEs.blastx \
-f 6 qseqid bitscore evalue stitle \
-F 15 -k 1 -e 1e-10 -p 8

# Obtain list of unique hits
awk -F '\t' '{print $1}' \
CarCar_VGPCheMyd_NoTEs.blastx \
| sort | uniq > CarCar_NoTEs_unique_hits.txt

echo $(wc -l CarCar_NoTEs_unique_hits.txt) # 121 unique hits

### Verify hits to non-TEs in CarCar dataset to check they are actual genes ### 

# Extract FASTA of the 121 CarCar seqids that hit CheMyd genes
module load samtools

samtools faidx CarCar-families.fa \
-r CarCar_NoTEs_unique_hits.txt > CarCar-gene-hits.fa

# Diamond blastx
$diamond blastx \
-d RepeatPeps.dmnd \
-q CarCar-gene-hits.fa \
-o CarCar_gene_hits_TE.1e10.blastx \
-f 6 qseqid bitscore evalue stitle \
-F 15 -k 1 -e 1e-10 -p 8

echo $(wc -l CarCar_gene_hits_TE.1e10.blastx) # 51 hits

# Make list of seqids that matched TEs 
awk -F '\t' '{print $1}' CarCar_gene_hits_TE.1e10.blastx \
> CarCar_gene_hits_TE_list.txt

# Reverse match TEs seqids -> obtain non-TEs only
grep -v -f CarCar_gene_hits_TE_list.txt \
CarCar_NoTEs_unique_hits.txt > CarCar_noTEs_unique_hits_checked.txt

# Check how many
echo $(wc -l CarCar_noTEs_unique_hits_checked.txt) # 70 hits

### Remove hits to non-TEs from CarCar-families.fa ###

# Make list of all seqids 
awk -F '\t' '{print $1}' \
CarCar-families.fa.fai > CarCarProteome_all_seqids.txt

echo $(wc -l CarCarProteome_all_seqids.txt) # 2261 lines

# Reverse match no TEs seqids -> obtain true TEs only
grep -v -f CarCar_noTEs_unique_hits_checked.txt \
CarCarProteome_all_seqids.txt > CarCar_true_TEs_seqids.txt

echo $(wc -l CarCar_true_TEs_seqids.txt) # 2191 lines

# Extract true TEs only from CarCar-families.fa (remove potential genes)
samtools faidx CarCar-families.fa  \
-r CarCar_true_TEs_seqids.txt > CarCar-families-filt.fa

# Check
samtools faidx CarCar-families-filt.fa
echo $(wc -l CarCar-families-filt.fa.fai) # 2191 entries
