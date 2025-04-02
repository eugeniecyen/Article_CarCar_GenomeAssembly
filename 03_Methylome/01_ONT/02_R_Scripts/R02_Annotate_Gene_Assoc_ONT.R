### R02_Annotate_Gene_Assoc.R ###

########################################################################################################################################

# Created by Charley, 2024

# Annotate and extract all CpG sites associated with genes for comparisons between
# ONT reference methylome and WGBS methylomes of 10 additional nesting females

########################################################################################################################################

###############################
###### Prep environment #######
###############################

####### Load packages ######

library(readxl)
library(tidyverse)
library(GenomicRanges)
library(genomation)
library(viridis)

####### Set directories ######

DIR <- file.path("/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes")
DIR_ONT <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/01_Working/02_ONT"
DIR_WGBS <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/01_Working/03_WGBS"
AnnoPath <- "/data/SBCS-EizaguirreLab/Turtle_Genome/10_Annotation/06_Functional_Annotation"

########################################################################################################################################

#################################################################
### Annotate and extract all sites associated with genes: ONT ###
#################################################################

####### Prep files ######

# Load Alice's BED file, already cleaned up from raw ONT file
# Chrom, position, coverage of CpG, modified C (hydroxy + meth, same as WGBS), fracMeth (% meth)
load(file.path(DIR_ONT, "myTurtleBedSum.RData")) 
head(myTurtleBedSum) # Check

# Load annotation
myannotGff3 <- rtracklayer::readGFF(file.path(AnnoPath, "gff3s/CarCar_QM_v1.21.12_Sc_Chr0_LongIso.gff3"))

# Load bed12 file
# NB Promoters are defined by options at genomation::readTranscriptFeatures function. 
# The default option is to take -1000,+1000bp around the TSS and you can change that. 
# following Heckwolf 2020 and Sagonas 2020, we consider 1500bp upstream and 500 bp downstream

# NB. ReadTranscriptFeatures defines the promoters. If you are using annotation with isoforms, you need to use unique.prom=FALSE option, otherwise
# promoter boundaries will not be assigned a gene name! 
# NB2. Even with longest isoforms only, CarCar annotation needs unique.prom to be set to FALSE.

myannotBed12 = readTranscriptFeatures(file.path(AnnoPath,"bed12s/CarCar_QM_v1.21.12_Sc_Chr0_LongIso.bed12"),
                                      remove.unusual = FALSE, up.flank = 1500, down.flank = 500, unique.prom=FALSE)

head(myannotBed12)

# Recursively change the gene names to keep only ID
getName <- function(x) {sub(";.*", "", sub(".*ID=", "", x))}

for (i in 1:length(myannotBed12)){
  myannotBed12[[i]]$name <- getName(myannotBed12[[i]]$name)
}

####### Convert ONT methylation bed file into a GRanges object #######

GRangeOBJ = makeGRangesFromDataFrame(
  data.frame(chr=myTurtleBedSum$chrom,
             start=myTurtleBedSum$start,
             end=myTurtleBedSum$start,
             fracMeth=myTurtleBedSum$fracMeth), keep.extra.columns = T)

####### Add feature type annotations to all sites ######

# Assign feature type to all sites first, so we can use it to guide finer-tuned annotation in the next step

# Annotate sites with genic parts
A = annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), feature = myannotBed12)

# Assign feature type to GRangeOBJ
# NB. Need to match names in myannotBed12, otherwise you get issues in next step when splitting annotation by genic vs intergenic
# names(myannotBed12) # Check

# Using info in A@members, replace rows where there are promoters, exons, introns or intergenic as appropriate
GRangeOBJ$featureType = ifelse(A@members[,1]==1, "promoters",
                               ifelse(A@members[,2]==1, "exons",
                                      ifelse(A@members[,3]==1, "introns", "intergenic")))

# Save: GRanges object with all feature type annotations, before filtering intergenic to <10kb away
saveRDS(GRangeOBJ, file = file.path(DIR_ONT, "GRanges_myTurtleBedSum_FeatureTypeAnnot_AllSites.RDS"))

####### Assign gene identities to all sites #######

## Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or, if intergenic, not further than 10 kb away from the TSS.

### Case 1: the feature is GENIC -> get annotation by intersection with bed file
# As site is on a gene (promoter, intron or exon), we now don't filter by distance from TSS for all DMS beforehand, otherwise 
# you could have a site at the end of a large gene that could be filtered out by accident

# Subset GRangesOBJ to retain DMS sites on genes only
GRangeOBJ1=GRangeOBJ[!GRangeOBJ$featureType %in% "intergenic"]

# Add empty column called geneInfo to fill with gene IDs
GRangeOBJ1$geneInfo <- NA

add_geneInfo_genic <- function(x, GRangeOBJ, annotBed12){
  ov = GenomicRanges::findOverlaps(
    annotBed12[[x]],
    GRangeOBJ[GRangeOBJ$featureType %in% x,])
  ## Add gene annotation to subject GRanges (i.e. left join)
  mcols(GRangeOBJ[GRangeOBJ$featureType %in% x,])[subjectHits(ov), "geneInfo"] = mcols(annotBed12[[x]])[queryHits(ov), "name"]
  return(GRangeOBJ)
}

myGRangeOBJ1 = add_geneInfo_genic("exons", GRangeOBJ1, myannotBed12) # Add gene names for sites on exons
myGRangeOBJ2 = add_geneInfo_genic("introns", myGRangeOBJ1, myannotBed12) # Add gene names for sites on introns
myGRangeOBJ_genic = add_geneInfo_genic("promoters", myGRangeOBJ2, myannotBed12) # Add gene names for sites on promoters

### Case 2: the feature is INTERGENIC: get the annotation by proximity to nearest TSS
# Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or,
# if intergenic, not further than 10 kb away from the TSS.

# Subset GRangesOBJ to retain intergenic sites only
GRangeOBJ2=GRangeOBJ[GRangeOBJ$featureType %in% "intergenic"]

a = annotateWithGeneParts(target = as(GRangeOBJ2,"GRanges"), feature = myannotBed12)

# Filter out sites that are further than 10 kb away from the TSS
rows2rm = which((a@dist.to.TSS$dist.to.feature>10000 | a@dist.to.TSS$dist.to.feature< -10000) &
                  rowSums(a@members) %in% 0)
if (is_empty(rows2rm)){  GRangeOBJ2 = GRangeOBJ2
} else { GRangeOBJ2 = GRangeOBJ2[-rows2rm,] }

# Re-annotate the subsetted object
b = annotateWithGeneParts(as(GRangeOBJ2,"GRanges"), myannotBed12)

# Get genes associated with these TSS
c = getAssociationWithTSS(b)

# Add these gene ID associations to GRangeOBJ2
GRangeOBJ2$geneInfo=c$feature.name

myGRangeOBJ_intergenic <- GRangeOBJ2

####### Merge the 2 cases back together into single GRangesOBJ and clean up #######

# Merge
myGRangeOBJ=c(myGRangeOBJ_genic, myGRangeOBJ_intergenic)

# Save
saveRDS(myGRangeOBJ, file = file.path(DIR_ONT, "GRanges_myTurtleBedSum_GeneInfoPerSite.RDS"))

# Clean up intermediate files
rm(GRangeOBJ1)
rm(GRangeOBJ2)
rm(myGRangeOBJ1)
rm(myGRangeOBJ2)
rm(myGRangeOBJ_genic)
rm(myGRangeOBJ_intergenic)
rm(a)
rm(A)
rm(b)
rm(c)
rm(GRangeOBJ)

########################################################################################################################################

##################################################################
### Annotate and extract all sites associated with genes: WGBS ###
##################################################################

####### Prep files ####### 

# Load methylKit uniteCov object
uniteCov <- readRDS(file.path(DIR_WGBS, "02_UniteCov_Objects/GenomePaper_WholeGenome_uniteCov75pc_Reloc3_Mothers.RDS")) 

# Create columns with metadata to include in final csv: names can't be those auto detected by makeGRangesFromDataFrame()
# These are included for checking the final functional annotation, due to previous issue where chrom assignment was mixed up
uniteCov$original_pos <- uniteCov$start
uniteCov$original_chrom <- uniteCov$chr

# Create vector of DMS with contig name and position in contig
# plus qvalue and meth diff to add as metadata -> can add column in final df
myPos = paste(uniteCov$chr, uniteCov$end, uniteCov$original_pos, uniteCov$original_chrom)

# Change the Pos vector into a dataframe): called df
df <- data.frame(chr=sapply(strsplit(myPos, " "), `[`, 1),
                 start=sapply(strsplit(myPos, " "), `[`, 2),
                 end=sapply(strsplit(myPos, " "), `[`, 2),
                 original_pos=sapply(strsplit(myPos, " "), `[`, 3),
                 original_chrom=sapply(strsplit(myPos, " "), `[`, 4))

# Convert into GRanges object
GRangeOBJ <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

####### Add feature type annotations to all sites ######

# Assign feature type to all sites first, so we can use it to guide finer-tuned annotation in the next step

# Annotate sites with genic parts
A = annotateWithGeneParts(target = as(GRangeOBJ,"GRanges"), feature = myannotBed12)

# Assign the feature type to GRangeOBJ
# NB. Need to match names in myannotBed12, otherwise you get issues in next step when splitting annotation by genic vs intergenic
# names(myannotBed12) # Check

# Using info in A@members, replace rows where there are promoters, exons, introns or intergenic as appropriate
GRangeOBJ$featureType = ifelse(A@members[,1]==1, "promoters",
                               ifelse(A@members[,2]==1, "exons",
                                      ifelse(A@members[,3]==1, "introns", "intergenic")))

# Save: GRanges object with all feature type annotations, before filtering intergenic to <10kb away
# This can be used for donut chart of where global CpG sites are distributed on genomic regions, where we want all intergenic regions
saveRDS(GRangeOBJ, file = file.path(DIR_WGBS, "03_GRanges_Objects/GRanges_75pc_Reloc3_Mothers_FeatureTypeAnnot_AllSites.RDS"))

####### Assign gene identities to all sites #######

## Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or, if intergenic, not further than 10 kb away from the TSS.

### Case 1: the feature is GENIC -> get the annotation by intersection with bed file
# As site is on a gene (promoter, intron or exon), we now don't filter by distance from TSS for all DMS beforehand, otherwise 
# you could have a site at the end of a large gene that could be filtered out by accident

# Subset GRangesOBJ to retain DMS sites on genes only
GRangeOBJ1=GRangeOBJ[!GRangeOBJ$featureType %in% "intergenic"]

# Add empty column to fill with gene IDs
GRangeOBJ1$feature.name <- NA

add_geneInfo_genic <- function(x, GRangeOBJ, annotBed12){
  ov = GenomicRanges::findOverlaps(
    annotBed12[[x]],
    GRangeOBJ[GRangeOBJ$featureType %in% x,])
  ## Add gene annotation to subject GRanges (i.e. left join)
  mcols(GRangeOBJ[GRangeOBJ$featureType %in% x,])[subjectHits(ov), "feature.name"] = mcols(annotBed12[[x]])[queryHits(ov), "name"]
  return(GRangeOBJ)
}

GRangeOBJ_ex = add_geneInfo_genic("exons", GRangeOBJ1, myannotBed12) # Add gene names for sites on exons
GRangeOBJ_ex_in = add_geneInfo_genic("introns", GRangeOBJ_ex, myannotBed12) # Add gene names for sites on introns
GRangeOBJ_genic = add_geneInfo_genic("promoters", GRangeOBJ_ex_in, myannotBed12) # Add gene names for sites on promoters

### Case 2: the feature is INTERGENIC: get the annotation by proximity to nearest TSS
# Heckwolf 2020: To be associated to a gene, the DMS had to be either inside the gene or,
# if intergenic, not further than 10 kb away from the TSS.

# Subset GRangesOBJ to retain intergenic sites only
GRangeOBJ2=GRangeOBJ[GRangeOBJ$featureType %in% "intergenic"]

a = annotateWithGeneParts(target = as(GRangeOBJ2,"GRanges"), feature = myannotBed12)

# Filter out sites that are further than 10 kb away from the TSS
rows2rm = which((a@dist.to.TSS$dist.to.feature>10000 | a@dist.to.TSS$dist.to.feature< -10000) &
                  rowSums(a@members) %in% 0)
if (is_empty(rows2rm)){  GRangeOBJ2 = GRangeOBJ2
} else { GRangeOBJ2 = GRangeOBJ2[-rows2rm,] }

# Re-annotate the subsetted object
b = annotateWithGeneParts(as(GRangeOBJ2,"GRanges"), myannotBed12)

# Get genes associated with these TSS
c = getAssociationWithTSS(b)

# Add these gene ID associations to GRangeOBJ2
GRangeOBJ2$feature.name=c$feature.name
GRangeOBJ_intergenic <- GRangeOBJ2


####### Merge the 2 cases back together into single GRangesOBJ and clean up #######

# Merge
myGRangeOBJ=c(GRangeOBJ_genic, GRangeOBJ_intergenic)

# Remove columns for original pos/chrom after checking it still matches!
myGRangeOBJ$original_pos <- NULL
myGRangeOBJ$original_chrom <- NULL

# Clean up intermediate files after checking all ok
rm(a)
rm(A)
rm(b)
rm(c)
rm(df)
rm(ov)
rm(rows2rm)
rm(myPos)
rm(list=ls(pattern="^GRangeOBJ"))

####### Subset uniteCov object for gene-associated sites only #######

# Subset uniteCov object for only sites associated to a gene (overlapping with a gene, <10kb from TSS for intergenic)
# Using methylKit::selectByOverlap with the GRanges object created above
uniteCov_Genes <- selectByOverlap(uniteCov, myGRangeOBJ)

# Save
saveRDS(uniteCov_Genes, file = file.path(DIR_WGBS, "03_GRanges_Objects/GenomePaper_uniteCov75pc_Reloc3_Mothers_GeneAssocSites.RDS"))


