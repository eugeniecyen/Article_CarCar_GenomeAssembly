#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=1G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley 2021, adapted from Chema
# Rename repeats with TEclass categories

#######################################################################

DIR=/data/scratch/btx902/01_Repeat_Masking/03_RepeatMasker/TEclass_renaming
cd $DIR

# Copy file over
cp ../../02_Filter_Bona_Fide_Genes/CarCar-families-filt.fa .

### Create list of missing entries ###

# Create tsv with columns: sequence name and TEclass
grep '>' CarCar-families-filt-TEclass.txt | \
sed 's/|/:/g' | sed 's/#/:/' | \
awk  ' BEGIN { FS = ":" } ; { OFS="\t" ; print $1,$4}' \
> table_of_changes.tsv

# Create list of sequence names in FASTA file
grep '>' CarCar-families-filt.fa | \
sed 's/|/:/g' | sed 's/#/:/' | \
awk  ' BEGIN { FS = ":" } ; { OFS="\t" ; print $1}' \
> fasta_entries

# Create list of sequence names in TEclass.txt file
grep '>' CarCar-families-filt-TEclass.txt | \
sed 's/|/:/g' | sed 's/#/:/' | \
awk  ' BEGIN { FS = ":" } ; { OFS="\t" ; print $1}' \
> text_entries

# Identify sequence names in FASTA file that are missing from TEclass.txt file
comm -13 <(sort text_entries) <(sort fasta_entries) > missing_entries

# Add missing entries to table_of_changes.tsv with classification "unclear"
cat missing_entries | awk '{ OFS="\t" ; print $1,"unclear"}' >> table_of_changes.tsv

### Replace existing TE classifications to those from TEclass in FASTA file ##

while read -r line
do
Target=$(echo $line | awk '{ print $1}')
rep=$(echo $line | awk '{ print $2}')
sed -i -E 's/'"$Target"'#.*/'"$Target"'#'"$rep"'/' CarCar-families-filt.fa
done < table_of_changes.tsv

# Rename file
mv CarCar-families-filt.fa  CarCar-families-filt-TEclass.fa

