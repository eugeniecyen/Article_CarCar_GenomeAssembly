### R01.1_PrepObjects_ChromArray_75pc_Reloc3_Adults.R ###

########################################################################################################################################

# Created by Charley, Jan 2024
# To be submitted on Apocrita with corresponding bash script as an array per chrommosome

# For 10 nesting females from Sal (Relocation 3, 2021 hatchery experiment)
# Create methylKit methylBase objects per chromosome, with sites covered in >75% of individuals

########################################################################################################################################

###############################
###### Prep environment #######
###############################

print(paste0("############### Started at ", Sys.time(), " ###############"))
start_time <- Sys.time()

####### Load packages ######

library(methylKit)
library(readxl)
library(plyr)
library(dplyr)
library(tidyverse)
library(GenomicRanges)


####### Set directories ######

dataPath <- "/data/SBCS-EizaguirreLab/Turtle_WGBS/04_Bismark_Methylation_Calls/destranded_methylation_calls"
metadataPath <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_WGBS_TSD_Genes/01_Working/01_Metadata"
OUTDIR_UniteCov <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_WGBS_TSD_Genes/01_Working/02_MethylKit_Objects/01_UniteCov_Objects/Mothers_Reloc3/by_chrom"

####### Pull arguments set in parent bash submission script ######

args <- commandArgs(trailingOnly = T)

# Set desired chromosome/s using task array ID, as specified in parent bash script
# Chrom number (i.e. array task ID) must be 1st arg in parent bash script
# Note this is SGE_TASK_ID - 1 because we have chromosomes 0-28 but the array must start at 1 -> (1:29)

print("Arg1 and Arg2:")
print(args)
chr = as.numeric(args[1])-1

# Set number of cores available (= 2nd argument ${NSLOTS} in parent bash script)
Num_cores = args[2]

print(paste0("##### Starting chromosome ", chr, " at ", Sys.time(), " #####")) 

########################################################################################################################################

###############################
# Create methylRawList object #
###############################

###### Add metadata ######

print("##### Reading metadata #####")

# Read in Excel file containing metadata
metadata <- readxl::read_xlsx(file.path(metadataPath, "Metadata_Mothers_Reloc3.xlsx"))

# methylKit only allows numeric treatments
metadata <- mutate(metadata,
                   trt_NUM = case_when(
                     Type == "Adult" ~ 1, 
                   ))

# Check
print("Printing head of metadata:")
head(metadata)
print("Printing no. of rows of metadata:")
nrow(metadata)

###### Load Bismark meth call files ######

# NB. .CpG_merged.cov.gz files were created following same pipeline here: 
# https://github.com/eugeniecyen/Article_CarCar_ThermalSublethalMeth/tree/master/Code/01_WGBS_Bismark_Pipeline 

print("##### Creating list of meth call files to analyse #####")

# Make a list of all destranded methylation call files
temp = list.files(path=dataPath,
                  pattern = paste0("chr", chr, ".CpG_merged.cov.gz" ),
                  full.names = T,
                  recursive = T)

# Subset list for only sample IDs that match those in metadata
temp <- temp[ grep(pattern=paste(metadata$Sample, collapse = "|"), x = temp) ]

# Check correct files are included in list
print("No. of files:")
length(temp) # check number of files

print("Printing head of file list:")
head(temp) # check 1st few files on the list

###### Make a methylRawList object ######

print("##### Creating methylRawList object with min coverage filter of 5X #####")

# Create object with min coverage filter of 8X
myobj.mincov8=methylKit::methRead(as.list(temp),
                                  pipeline='bismarkCoverage',
                                  mincov=8,
                                  sample.id=as.list(metadata$Sample),
                                  assembly="CarCar_QM_v1_2021_12_Scaff_With_Chr0.fasta",
                                  treatment=metadata$trt_NUM,
                                  context="CpG")

print("##### Finished #####")

########################################################################################################################################

###############################
## Filtering and normalising ##
###############################

####### Filter based on coverage ######

print("##### Filtering methylRawList object by maximum >99.9% and minimum <5X coverage #####")

filt.myobj.mincov8=filterByCoverage(myobj.mincov8, lo.count=8, lo.perc=NULL,
                                    hi.count=NULL, hi.perc=99.9)

print("##### Finished #####")

####### Normalise coverage #######

# Not necessarily needed if coverage is similar across samples, but good to do in case
# Together with removing extreme coverage, will help reduce the bias in  statistical tests 
# that might occur due to systematic over-sampling of reads in certain samples

print("##### Normalising coverage #####")

normFil.myobj.mincov8=normalizeCoverage(filt.myobj.mincov5)

print("##### Finished #####")

####### Remove files that are no longer needed to save memory #######

rm(myobj.mincov8)
rm(filt.myobj.mincov8)

########################################################################################################################################

#####################
### Merge samples ###
#####################

#######  Create methylBase object ####### 

print("##### Creating uniteCovALL object with CpG present in >75% of individuals per treatment #####")

uniteCov75pc = methylKit::unite(normFil.myobj.mincov5, destrand=FALSE, mc.cores = Num_cores, 
                                min.per.group=as.integer(length(metadata$Sample)*0.75))

uniteCov75pc = as(uniteCov75pc,"methylBase")

# Check
print("Printing head of uniteCov75pc:")
head(uniteCov75pc)
print("No. of rows in uniteCov75pc is:")
nrow(uniteCov75pc)

# Save
print("##### Saving #####")

saveRDS(uniteCov75pc, file = file.path(OUTDIR_UniteCov, paste0("chr", chr, "_uniteCov75pc.RDS")))

print("##### Finished saving #####")

########################################################################################################################################

end_time <- Sys.time()

print(paste0("############### All finished for chromosome ", chr, " at ", end_time, " ###############"))

runtime <- end_time - start_time

print("Total run time:")
print(runtime)














