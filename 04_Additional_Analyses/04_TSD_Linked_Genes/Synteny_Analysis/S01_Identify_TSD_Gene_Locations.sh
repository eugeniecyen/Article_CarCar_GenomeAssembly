#!/bin/bash
#$ -cwd
#$ -pe smp 1      
#$ -l h_rt=1:0:0
#$ -l h_vmem=4G
#$ -j y

#######################################################################

# Created by: James Gilbert (https://github.com/evoJDG), 2024

###### Identify Locations of SD Genes in loggerhead (Ccar) and hawksbill (Eimb) genomes
# We want to know whether genes related to sex determination are on the same chromosomes
# in loggerhead, green, leatherback and hawksbill turtles

# 1. Download genome for green turtle
# 2. Extract sequences for SD genes (coordinates in metadata)
# 3. Loop for Ccar and Eimb:
# 		BLAST sequences against genome
#		Verify genes are correct (look up names)
#		Extract locations of SD genes and add to table of locations

#######################################################################

############################
### Prepare Environment ####
############################

## load modules
module load samtools/1.10
module load blast+/2.11.0

module load anaconda3/2020.02
conda activate Orthofinder

## Set directories
# for the project
DIR=/data/SBCS-EizaguirreLab/James/SD_Flex/Submission-2
SYNDIR=${DIR}/01_Synteny
METADIR=${DIR}/Genome_SD_Flex/00_Metadata
SCRIPTDIR=${DIR}/Genome_SD_Flex/01_Synteny

# for assembly genome and annotation
ASSEMDIR=/data/SBCS-EizaguirreLab/Turtle_Genome
ANNOTDIR=${ASSEMDIR}/10_Annotation/06_Functional_Annotation
GENOMEDIR=${ASSEMDIR}/00_Final_Assemblies
## metadata
# The locations of the SD genes were published in Bentley et al. metadata
META=${METADIR}/Bentley_TSDgene_info.txt

########################################
### 1. Download and Prepare Genomes ####
########################################

### Paths for CC, Download genomes for CM and DC
# and find longest transcripts for Orthofinder
# will save the green genome here:
CM_DIR=${SYNDIR}/data/Genomes/Cmydas
CC_DIR=${SYNDIR}/data/Genomes/Ccaretta
DC_DIR=${SYNDIR}/data/Genomes/Dcoriacea
EI_DIR=${SYNDIR}/data/Genomes/Eimbricata
mkdir -p $CM_DIR $CC_DIR/blast $DC_DIR $EI_DIR

##### Loggerhead #####
# WG fasta
CC_REF=${GENOMEDIR}/CarCar_QM_v1_2021_12_Scaff_With_Chr0.fasta 

# Extract CDS
CC_REF_CDS=${CC_DIR}/CarCar_QM_v1_2021_12_Scaff_With_Chr0.cds.fa

#### Annotations
CC_GFF=${ANNOTDIR}/gff3s/CarCar_QM_v1_2021_12_Scaff_Func_Anno.gff3

### the Cc gff doesn't match genome. Filter gff for only seqs that map to full chromosomes i.e. those with SLK_ragtag_063 prefix.
CC_GFF_filtered=${CC_DIR}/CarCar_QM_v1_2021_12_Scaff_Func_Anno.filtered.gff3
awk '$1 ~ /^SLK063_ragtag/ || $1 ~ /^#/' $CC_GFF > $CC_GFF_filtered

# use gffread to extract CDS
gffread -x ${CC_REF_CDS} \
-g ${CC_REF} \
${CC_GFF_filtered}

#### Longest isoform annotation
CC_GFF_longestIso=${ANNOTDIR}/gff3s/CarCar_QM_v1_2021_12_Scaff_Func_Anno_LongestIso.gff3

### the Cc gff doesn't match genome. Filter gff
CC_GFF_longestIso_filtered=${CC_DIR}/CarCar_QM_v1_2021_12_Scaff_Func_Anno.LongestIso.filtered.gff3

awk '$1 ~ /^SLK063_ragtag/ || $1 ~ /^#/' $CC_GFF_longestIso > $CC_GFF_longestIso_filtered

##### Extract longest isoforms of CDS
# output for long Iso of CDS
CC_REF_CDS_longestIso=${CC_DIR}/CarCar_QM_v1_2021_12_Scaff_With_Chr0.LongestIso.cds.fa

### use gffread to extract CDS
gffread -x ${CC_REF_CDS_longestIso} \
-g ${CC_REF} \
$CC_GFF_longestIso_filtered

##### Green #####

echo -e "\n\n\n##### downloading C. mydas genome #####\n\n\n"

# the assembly is GCF_015237465.2, available at the URL below

# shortcut to reference assembly. Note we are using this version, not longest tx, to be consistent with the gene locations in Bentley et al metadata (SI7)
CM_REF=${CM_DIR}/ncbi_dataset/data/GCF_015237465.2/GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna

cd $CM_DIR

# download the genome from ncbi
curl -OJX GET https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_015237465.2/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT \
-H "Accept: application/zip"

# unzip and remove the download
unzip ncbi_dataset.zip 
rm ncbi_dataset.zip

# find longest transcripts
# script is downloaded from Orthofinder
python ${SCRIPTDIR}/primary_transcript.py ${CM_DIR}/ncbi_dataset/data/GCF_015237465.2/cds_from_genomic.fna

##### Leatherback #####
echo -e "\n\n\n##### downloading D. coriacea genome #####\n\n\n"

# the assembly is GCF_015237465.2, available at the URL below
cd $DC_DIR

# download the genome from ncbi
curl -OJX GET https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_009764565.3/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT \
-H "Accept: application/zip"

# unzip and remove the download
unzip ncbi_dataset.zip 
rm ncbi_dataset.zip 

# find longest transcripts
# script is downloaded from Orthofinder
python ${SCRIPTDIR}/primary_transcript.py ${DC_DIR}/ncbi_dataset/data/GCF_009764565.3/cds_from_genomic.fna

#### hawksbill #####

echo -e "\n\n\n##### downloading E. imbricata genome #####\n\n\n"

# Hawksbill downloaded from figshare, not ncbi

cd $EI_DIR

### cds and annotation here
wget -O ${EI_DIR}/downloaded_file.tar.gz "https://figshare.com/ndownloader/files/41765319"

tar -xvzf ${EI_DIR}/downloaded_file.tar.gz
rm ${EI_DIR}/downloaded_file.tar.gz

# find longest transcripts
# script is downloaded from Orthofinder
python ${SCRIPTDIR}/primary_transcript.py ${EI_DIR}/NSCA_data/Eretmochelys_imbricata.Genome.V1.cds

#### Annotation
EI_GFF=${EI_DIR}/NSCA_data/Eretmochelys_imbricata.Genome.V1.gff3

########################################
## 2. Extract Sequences for SD Genes ###
########################################

echo -e "\n\n\n##### extracting SD sequences #####\n\n\n"

# file to save the SD sequences to. create empty file to append seqs to
SDSEQS=${CM_DIR}/ncbi_dataset/data/GCF_015237465.2/SD_Sequences.fna
> $SDSEQS

# read line was being awkward so looping through an index
NROW=$( sed -n '$=' $META ) # to loop through
MAP=$( tail -n +2 ${METADIR}/CM_Chrom_Map.tsv ) # mapping chromosome names from Bentley SI to genome files

# loop through the lines
for i in $(eval echo "{2..$NROW}"); do

	GENE=$( cut -f 1 $META | sed -n ${i}p )
	CHRNUM=$( cut -f 11 $META | sed -n ${i}p | sed 's/SUPER_//g' ) # just want the chromosome number
	
	# In Bentley metadata, chromosomes are named SUPER_1, SUPER_2 etc, but in the sequence file, they have sequence names e.g. NC_057849.1
	# $MAP is the file mapping chromosome number (1,2,3) to sequence ID (NC_12345.1)
	# use the chromosome number to extract the corresponding sequence name
	CHR=$( echo "$MAP" | sed -n ${CHRNUM}p | cut -f 10 )
	START=$(  cut -f 12  $META | sed -n ${i}p | sed 's/[",]//g' )
	END=$(  cut -f 13  $META | sed -n ${i}p | sed 's/[",]//g' )

	# cut the corresponding sequence out of the genome and print to the file to be blasted
	samtools faidx $CM_REF ${CHR}:${START}-${END} | sed "/>/c\>$GENE" \
	| awk '/^>/ { if (NR > 1) print ""; printf("%s\n",$0); next } { printf("%s",$0) } END { printf("\n") }' \
	>> $SDSEQS

done

##### Find SD genes in CC and EI

### Will loop through - need the following info in each loop
Species=("CC" "EI")
CDS_in=("${CC_DIR}/CarCar_QM_v1_2021_12_Scaff_With_Chr0.cds.fa" "${EI_DIR}/NSCA_data/Eretmochelys_imbricata.Genome.V1.cds" )
outDir=("${CC_DIR}" "${EI_DIR}")
GFF=("$CC_GFF" "$EI_GFF")

# Note Spec iterable - don't want to clash with any $i in the scripts called from the loop
for Spec in {0..1}; do

	echo ${Species[$Spec]}
    echo ${CDS_in[$Spec]}
    echo ${outDir[$Spec]}

    ############################################
    ###  3. BLAST sequences against genome #####
    ############################################

    ### Here we BLAST the SD gene sequences against the CC and EI genomes
    # ### We do it in a loop using the S02_blastSDseqs.sh

    source ${SCRIPTDIR}/utils/blastSDseqs.sh

    source ${SCRIPTDIR}/utils/addGeneLocationToFile.sh

done

##### Compare results from searches in CC and EI genomes. Make sure it looks roughly even
cut -f7 ${SYNDIR}/data/Locations/SDGene_Locations_CC.txt | sort | uniq -c

cut -f7 ${SYNDIR}/data/Locations/SDGene_Locations_EI.txt | sort | uniq -c