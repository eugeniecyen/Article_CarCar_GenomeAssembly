#!/bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=4:0:0
#$ -l h_vmem=4G
#$ -j y

#######################################################################

# Created by: James Gilbert (https://github.com/evoJDG), 2024

###### Run Orthofinder to find orthologous genes across the species in the study
#### We have the coding sequences for each species' genome
#### And use Orthofinder to find single-copy homologues present in all species
#### For each orthogroup, Orthofinder produces an unaligned multi-species fasta

#### Note Orthofinder can do various other steps, but we stop at this point
#### using the -os flag ("Stop after writing sequence files for orthogroups (requires '-M msa')")

###### After identifying orthologues, we want to make sure none of them overlap with SD genes
#### Blast the orthologues against all SD seqs (from green turtle) and exclude any orthologues with hits

##### Finally, we look up the locations of the OGs in the Car car genome for downstream analyses

#######################################################################

#################################
######## Prep Environment  ######
#################################

#### load anaconda and activate orthofinder env (where orthofinder is installed)
module load anaconda3/2020.02 
conda activate Orthofinder # v2.5.4 of orthofinder is installed

## set directories
DIR=/data/SBCS-EizaguirreLab/James/SD_Flex/Submission-2
SYNDIR=${DIR}/01_Synteny
METADIR=${DIR}/Genome_SD_Flex/00_Metadata
SCRIPTDIR=${DIR}/Genome_SD_Flex/01_Synteny/
ORTHDIR=${DIR}/02_Methylation

ORTHRESULTSDIR=${ORTHDIR}/OrthoFinder/Results_Dec14 ###### IMPORTANT: THIS MAY NEED TO BE UPDATED FOR LATER ORTHOFINDER RUNS ######
mkdir -p $ORTHDIR/blast

# for assembly genome and annotation
ASSEMDIR=/data/SBCS-EizaguirreLab/Turtle_Genome
ANNOTDIR=${ASSEMDIR}/10_Annotation/06_Functional_Annotation
GENOMEDIR=${ASSEMDIR}/00_Final_Assemblies

CM_DIR=${SYNDIR}/data/Genomes/Cmydas

# for getting gene locations
CC_GFF=${ANNOTDIR}/gff3s/CarCar_QM_v1_2021_12_Scaff_Func_Anno.gff3

##### IMPORTANT ##### Don't re run Orthofinder every time - it will write to a new dir and the script won't work
# Control by setting this parameter to Y (yes) or N (no):
RUNORTH=N

if [[ "$RUNORTH" == "Y" ]]; then 

	#### Copy over the longest transcript genomes
	cp ${SYNDIR}/data/Genomes/Ccaretta/CarCar_QM_v1_2021_12_Scaff_With_Chr0.LongestIso.cds.fa \
	${ORTHDIR}/Ccaretta.cds_from_genomic_longestIsoform.fa

	cp ${SYNDIR}/data/Genomes/Cmydas/ncbi_dataset/data/GCF_015237465.2/primary_transcripts/cds_from_genomic.fna \
	${ORTHDIR}/Cmydas.cds_from_genomic_longestIsoform.fa

	cp ${SYNDIR}/data/Genomes/Dcoriacea/ncbi_dataset/data/GCF_009764565.3/primary_transcripts/cds_from_genomic.fna \
	${ORTHDIR}/Dcoriacea.cds_from_genomic_longestIsoform.fa

	cp ${SYNDIR}/data/Genomes/Eimbricata/NSCA_data/primary_transcripts/Eretmochelys_imbricata.Genome.V1.cds \
	${ORTHDIR}/Eimbricata.cds_from_genomic_longestIsoform.fa

	#################################
	######### Run Orthofinder  ######
	#################################

	echo "starting Orthofinder"

	# -t ${NSLOTS} uses available threads
	# running orthofinder
	# -S diamond option is the sequence search option
	# diamond is faster than blast
	# -os = Stop after writing sequence files for orthogroups. we don't need alignments, just sequences and names.
	# -m msa is for inferring trees and is required for -os
	# -d indicates DNA (nucleotide data)
	orthofinder -S diamond -t ${NSLOTS} -M msa -os -f ${ORTHDIR} -d

	echo -e "**** FINISHED ORTHOFINDER *****"

fi # finish if run orthofinder

#################################
#### Look for OGs in SD genes ###
#################################

#### IMPORTANT ##### Don't re run this every time - it takes ages
# Control by setting this parameter to Y (yes) or N (no):
Find_OGs_in_SD=Y

# set files to store output (outside of if so can be used further down, even if not looking for OGs in SD genes)
OG_USE=${ORTHRESULTSDIR}/blast/OGs_forAnalysis.txt
OG_EXCLUDE=${ORTHRESULTSDIR}/blast/OGs_toExclude.txt

### We loop through the OG sequences created by Orthofinder
### And blast against all the SD seqs for green turtle provided in Bentley's metadata. 
### We only assign the OG to the OGs_forAnalysis.txt file if there are not hits with the SD genes.
### if there is a hit against any SD seq, we assign to the OGs_toExclude list.

### Finally, we want the locations of the orthologous genes in the Car car genome.

if [[ "$Find_OGs_in_SD" == "Y" ]]; then

	echo -e "***** LOOKING FOR SD GENES IN ORTHOGROUPS *****"

	# make empty files to store output
	> $OG_USE
	> $OG_EXCLUDE

	## Identify OGs that are also TSD genes by blasting each OG.fa against a file containing the seqs
	# for all species for all TSD genes:
	SDSEQS=${CM_DIR}/ncbi_dataset/data/GCF_015237465.2/SD_Sequences.fna

	# the list of single copy orthogroups is all the files in here 
	OGS=$( ls ${ORTHRESULTSDIR}/Single_Copy_Orthologue_Sequences | sed 's/\.fa//g' )

	# loop through OG fastas
	for OG in $OGS; do

		# $OG is the name of orthogroup, this is full file path:
		OGfile=${ORTHRESULTSDIR}/Single_Copy_Orthologue_Sequences/${OG}.fa

		# output file for blast query
		mkdir -p ${ORTHRESULTSDIR}/blast
		OUT_OGBLAST=${ORTHRESULTSDIR}/blast/$OG.blast.out

		# blastn (nt vs nt) the query OG sequence
		# against the file containing seqs for all genes and species
		# save all outputs to a dir created in the Ortho results
		blastn -query $OGfile \
		-subject $SDSEQS \
		-outfmt 6 -out $OUT_OGBLAST

		# print OG name to either list of OGs to use or exclude

		# if the blast output has no size (! -s) i.e. there are not matches in TSD genes
		# print to list for use in downstream steps
		if [[ ! -s ${OUT_OGBLAST} ]]; then echo $OG >> $OG_USE; fi

		# if the blast output is not empty (has size; -s) i.e. there are matches in TSD genes
		# print to list for exclusion from downstream steps
		# also add which gene the OG corresponds to
		if [[ -s ${OUT_OGBLAST} ]]; then echo $OG >> $OG_EXCLUDE; fi

	done # end loop to find OGs that are TSD genes



	#### Get locations for Orthogroups ####

	echo -e "***** FINDING OG LOCATIONS IN THE LOGGERHEAD GENOME *****"

	# make empty results file
	ORTHLOCS=${SYNDIR}/data/Locations/Orthogroup_Locations.txt
	echo -e "gene\tsequence\tchromosome\tstart\tend" > $ORTHLOCS

	## this is where the SD genes are stored - we do a final check to make sure no OGs are also in the SD file
	SDLOCS=${SYNDIR}/data/Locations/SDGene_Locations_CC.txt

	if [[ ! -f "$SDLOCS" ]]; then 
		echo "Error: File '$SDLOCS' does not exist." >&2
		exit 1	
	fi

	### here we look up the locations of OG sequences in the CC genome.
	# loop through OG fastas
	while read OG; do

		# the sequences for each OG is in a fasta file:
		OGfile=${ORTHRESULTSDIR}/Single_Copy_Orthologue_Sequences/${OG}.fa

		# find the name of the CC sequence which we will use to search the SD sequences
		SEQ=$( grep "CarCar" $OGfile | sed 's/>//g' )

		### important thing: if the sequence is an SD gene that has not been picked up by the blast
		# this is the case for PDGFA, we don't want its locations here
		# grep q will return true or false
		if grep -q $SEQ $SDLOCS ; then

			# if the sequence is in the SD locations
			# remove it from the OG_USE list
			# and add to the OG exclude list
			sed -i "/$OG/d" $OG_USE

			echo $OG >> $OG_EXCLUDE

		else

			### look up the sequence in the annotation
			GENEINFO=$(grep $SEQ $CC_GFF)

			# extract chromosome, start and end from respective columns from the annotation file
			CHRM=$( echo "$GENEINFO" | cut -f 1 | head -n 1 )
			START=$( echo "$GENEINFO" | cut -f 4 | sort | head -n 1 )
			END=$( echo "$GENEINFO" | cut -f 5 | sort | tail -n 1 )

			# add the details to the output file
			echo -e "$OG\t$SEQ\t$CHRM\t$START\t$END" >> $ORTHLOCS

		fi


	done < "$OG_USE"

fi # finish if check OGs in SD

