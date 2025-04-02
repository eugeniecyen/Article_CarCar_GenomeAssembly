### R01.2_MergeChrom_75pc_Reloc3_Adults.R ###

########################################################################################################################################

# Created by Charley, Jan 2024
# To be submitted on Apocrita with corresponding bash script

# For 10 nesting females from Sal (Relocation 3, 2021 hatchery experiment)
# Merge together uniteCov75pc files by chromosome (produced in array by chromosome script) into single, whole genome file

########################################################################################################################################

###############################
###### Prep environment #######
###############################

print(paste0("############### Started at ", Sys.time(), " ###############"))
start_time <- Sys.time()

####### Pull arguments set in parent bash submission script ######

args <- commandArgs(trailingOnly = T)

# Set number of cores available (${NSLOTS} in parent bash script)
Num_cores = args[1]
print("Num_cores:")
print(args)

####### Load packages ######

library(methylKit)
library(readxl)
library(plyr)
library(dplyr)
library(tidyverse)
library(GenomicRanges)

####### Set directories ######

dataPath <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_WGBS_TSD_Genes/01_Working/02_MethylKit_Objects/01_UniteCov_Objects/Mothers_Reloc3/by_chrom"
OUTDIR_UniteCov <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_WGBS_TSD_Genes/01_Working/02_MethylKit_Objects/01_UniteCov_Objects/Mothers_Reloc3"

########################################################################################################################################

################################
# Merge chrom uniteCov objects #
################################

### Load in chromosome uniteCov files ###

print("##### Combining chromosome uniteCov75pc_WholeGenome files into one dataframe #####")

List_uniteCov <- list()
List_dfs_uniteCov <- list()

# Run loop adding uniteCov file for each chromosome into List_uniteCov
# NB. R cannot handle starting loop on 0 -> need to bypass
# Create list of chrom, with 0 at end to match order in ref genome - solves object issues later: has to be in same order as ref genome
chrom_number <- 1:28
chrom_number <- append(chrom_number, '0', after=28)
chroms <- c(paste0("chr", chrom_number))

for (chr in chroms) {
  print(chr)
  List_uniteCov[[ chr ]] <- readRDS(file.path(dataPath, paste0(chr, "_uniteCov75pc.RDS")))
  List_dfs_uniteCov[[ chr ]] <- getData(List_uniteCov[[ chr ]])
}

# Bind together into a data frame
df_uniteCov_AllChrms <- do.call(rbind, List_dfs_uniteCov)

print("##### Finished #####")

### Convert dataframe into a methylBase object ###

# Pull over metadata to fill the "slots" of the methylBase object
print("##### Creating methylKit object uniteCov75pc_WholeGenome from dataframe #####")

# Added line adding coverage.index to JG's code = missing metadata needed for methylKit::reorganize() to work
uniteCov_WholeGenome <- new(Class = "methylBase", rbind(df_uniteCov_AllChrms),
                            sample.ids=List_uniteCov[[ 1 ]]@sample.ids,
                            destranded=List_uniteCov[[ 1 ]]@destranded,
                            assembly=List_uniteCov[[ 1 ]]@assembly,
                            context=List_uniteCov[[ 1 ]]@context,
                            treatment=List_uniteCov[[ 1 ]]@treatment,
                            resolution=List_uniteCov[[ 1 ]]@resolution,
                            numCs.index=grep("numCs", colnames(df_uniteCov_AllChrms)), # need to tell methylkit which columns numCs
                            numTs.index=grep("numTs", colnames(df_uniteCov_AllChrms)), # and which columns contain numTs
                            coverage.index=grep("coverage", colnames(df_uniteCov_AllChrms)) # and which columns contain coverage
)

print("##### Finished #####")

print("Printing head of uniteCov75pc_WholeGenome:")
head(uniteCov_WholeGenome)
print("No. of rows in uniteCov75pc_WholeGenome is:")
nrow(uniteCov_WholeGenome)

### Save uniteCov_WholeGenome file ###

print("##### Saving #####")

saveRDS(uniteCov_WholeGenome, file.path(OUTDIR_UniteCov, "GenomePaper_WholeGenome_uniteCov75pc_Reloc3_Mothers.RDS"))

print("##### Finished #####")

########################################################################################################################################

print(paste0("############### All finished at ", Sys.time(), " ###############"))

end_time <- Sys.time()
runtime <- end_time - start_time

print("Total run time:")
print(runtime)
