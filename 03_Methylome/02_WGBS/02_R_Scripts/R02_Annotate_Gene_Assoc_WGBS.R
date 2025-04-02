### R02_Annotate_Gene_Assoc_WGBS.R ###

########################################################################################################################################

# Created by Charley, 2024
# Annotate and extract all CpG sites associated with genes in WGBS methylomes (sites present in >75% of individuals)

########################################################################################################################################

###############################
###### Prep environment #######
###############################

####### Load packages ######

library(methylKit)
library(readxl)
library(tidyverse)
library(GenomicRanges)
library(genomation)

####### Set directories ######

DIR_WGBS <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/01_Working/03_WGBS"
AnnoPath <- "/data/SBCS-EizaguirreLab/Turtle_Genome/10_Annotation/06_Functional_Annotation"

########################################################################################################################################

##################################################################
### Annotate and extract all sites associated with genes: WGBS ###
##################################################################

####### Prep files ####### 

# Load uniteCov object
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

# Subset uniteCov object for gene-associated sites
# Using methylKit::selectByOverlap with the GRanges object created above
uniteCov_Genes <- selectByOverlap(uniteCov, myGRangeOBJ)

# Save
saveRDS(uniteCov_Genes, file = file.path(DIR_WGBS, "03_GRanges_Objects/GenomePaper_uniteCov75pc_Reloc3_Mothers_GeneAssocSites.RDS"))
