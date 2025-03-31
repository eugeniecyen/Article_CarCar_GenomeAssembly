### R01_methylStat_methylKit.R ###

########################################################################################################################################

# Created by Alice Balard, 2023 (updated in 2025)
# Generate ONT methylation call stats and merge 5mC/5hmC sites, since WGBS cannot distinguish between them

########################################################################################################################################

###############################
###### Prep environment #######
###############################

####### Load packages ######

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
##BiocManager::install("methylKit")
##BiocManager::install("genomation")
##BiocManager::install("plyranges")

library(methylKit)
library(tidyr)
library(dplyr)
library(data.table)
library(genomation)
library(rlang)
library(plyranges)
library(ggplot2)
library(ggrepel)

########################################################################################################################################

#################################
###### Generate meth file #######
#################################

rerun = FALSE
if (rerun == FALSE){
    load(file = "/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/myTurtleBedSum_2025.RData")
  paste("Number of modified Cs:")
  nrow(myTurtleBedSum)
  } else {
    ## 1. Reading the methylation calls from sorted Guppy alignments (should be bam/sam like Bismark)
    ## Use fread, faster than read.csv (NB: tab delimited in bash)
    myTurtleBed = fread("/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/bambeds/turtleONTmerged.sorted.bam_combineStrands_tabdelim.bed", 
                        sep = "\t", header = F)
    ## Clean bed file
    ## https://github.com/nanoporetech/modkit
    names(myTurtleBed) = c("chrom", "start", "end", "modbase", "score",
                           "strand", "start2", "end2", "color",
                           "Nvalid_cov", "fracMod", "Nmod", "Ncanonical","Nother_mod",
                           "Ndelete", "Nfail", "Ndiff", "Nnocall")
    
    ##########################
    ## A. Descriptive stats ##
    ##########################
    ## Calculate methylation statistics
    
    ## number of 5hmc
    # Npos_5hmC = nrow(myTurtleBed[myTurtleBed$modbas == "h" & myTurtleBed$Nmod > 0,])
    # ## number of 5mc
    # Npos_5mC = nrow(myTurtleBed[myTurtleBed$modbas == "m" & myTurtleBed$Nmod > 0,])
    
    ## number of positions with 5hmc ONLY (=no pair reads with 5mc)
    Npos_5hmC_ONLY = nrow(myTurtleBed[myTurtleBed$modbas == "h" & myTurtleBed$Nmod > 0 & myTurtleBed$Nother_mod == 0,])
    Npos_5hmC_ONLY ## 120,983
    
    ## number of positions with 5mc ONLY (=no pair reads with 5hmc)
    Npos_5mC_ONLY = nrow(myTurtleBed[myTurtleBed$modbas == "m" & myTurtleBed$Nmod > 0 & myTurtleBed$Nother_mod == 0,])
    Npos_5mC_ONLY ## 22,327,230
    
    ## number of positions with 5mc AND 5hmc
    Npos_mCANDhmC = nrow(myTurtleBed[myTurtleBed$modbas == "m" & myTurtleBed$Nmod > 0 & myTurtleBed$Nother_mod > 0,])
    Npos_mCANDhmC ## 2,986,762
    
    ## number of positions with only unmethylated C
    Npos_nonmethC_only = nrow(myTurtleBed[myTurtleBed$modbas == "m" & myTurtleBed$Nmod == 0 & myTurtleBed$Nother_mod == 0,])
    Npos_nonmethC_only # 1,014,100
    
    ## Sanity check
    Npos_5hmC_ONLY + Npos_5mC_ONLY + Npos_mCANDhmC + Npos_nonmethC_only # 26,449,075
    
    ## total number of C sequenced
    nrow(myTurtleBed)/2 # 26,449,075
    
    ## Simplify format -> merge h and m
    myTurtleBedSum = myTurtleBed[ ,.(coverage = sum(Nvalid_cov), modC = sum(Nmod, Nother_mod)), by = .(chrom, start, end)]
    myTurtleBedSum$fracMeth = myTurtleBedSum$modC/myTurtleBedSum$coverage
    
    ## Save RData
    save(myTurtleBedSum, file = "/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/myTurtleBedSum.RData")
}

## Percentage of methylated (>0.7) sites
nrow(myTurtleBedSum[myTurtleBedSum$fracMeth >= 0.7,])/nrow(myTurtleBedSum)
## 74.13% of the sites are methylated above 0.7

## Percentage of methylated (>0.99) sites
# nrow(myTurtleBedSum[myTurtleBedSum$fracMeth >= 0.99,])/nrow(myTurtleBedSum)
# ## 18.71% of the sites are methylated above 0.99
