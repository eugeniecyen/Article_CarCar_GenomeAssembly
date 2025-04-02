### R02.2_Check_TSD_Gene_Locations_EI.R ###

########################################################################################################################################

# Created by: James Gilbert (https://github.com/evoJDG), 2024

# S01, S02 and S03 have found the positions of the SD-related genes in CC and EI, and mapped them to the same genes in CM and DC genomes.
# Here, we plot the positions to show the synteny and highlight non-syntenic cases.
# We plot only for CC against each of the other 2.

########################################################################################################################################

###########################
##### Prep Environment ####
###########################

library(dplyr)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library( viridis )
library(viridisLite)
library( gtools )

# Load data
DIR <- "/data/SBCS-EizaguirreLab/James/SD_Flex/Submission-2"
Plot_DIR = file.path(DIR, "01_Synteny/Plots"); dir.create(Plot_DIR, showWarnings = F)

Locations_All_Plot = read.table( file.path(DIR, "01_Synteny/data/Locations/SDGene_Locations_All_Plot.txt"), sep = "\t",
                                 header = T)

########################################################################################################################################

#########################
#### Plot Locations #####
#########################

##### Prepare data #####

### chromosome data ###

# need chromosome name, start coord (always 0) and length
# note the order of the chromosomes in the dataframes and angle of the plot determine where they are located
# the current orientation splits the circle down the middle
# with chromosome 1 at the top and 28 at the bottom. CC on the left and others on the right

# chrom data for CC
CC_chrom_data <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/CC_Chrom_Map.tsv"), sep = "\t", header = F ) %>%
  select( V1, V2) %>%
  rename( Chromosome.name = V1, Seq.length = V2) %>%
  mutate(Seq.length = as.numeric(Seq.length), # must be numeric
         Chromosome.name = gsub("SLK063_ragtag_chr", "CC_", Chromosome.name),
         start = 0 ) %>%
  select( Chromosome.name, start, Seq.length ) %>%
  filter( Chromosome.name != "CC_0" ) %>%
  group_by(Chromosome.name) %>% 
  arrange(rev(.)) %>% as.data.frame

# for CM
CM_chrom_data <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/CM_Chrom_Map.tsv"), sep = "\t", header = T) %>%
  filter(grepl("NC_", RefSeq.seq.accession),
         Chromosome.name %in% 1:28)   %>%
  mutate( Chromosome.name = paste0("CM_", Chromosome.name),
          start = 0) %>%
  select(Chromosome.name, start,  Seq.length)

# for DC
DC_chrom_data <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/DC_Chrom_Map.tsv"), sep = "\t", header = T) %>%
  filter(grepl("NC_", RefSeq.seq.accession),
         Chromosome.name %in% 1:28)   %>%
  mutate( Chromosome.name = paste0("DC_", Chromosome.name),
          start = 0) %>%
  select(Chromosome.name, start,  Seq.length)

# for EI
EI_chrom_data <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/EI_Chrom_Map.tsv"), sep = "\t", header = T) %>%
  filter(grepl("^CM", GenBank.seq.accession),
         Chromosome.name %in% 1:28)   %>%
  mutate( Chromosome.name = paste0("EI_", Chromosome.name),
          start = 0) %>%
  select(Chromosome.name, start,  Seq.length)

# Combine to single df for each pair 
CC_CM_chrom_data <- rbind( CM_chrom_data, CC_chrom_data)
CC_DC_chrom_data <- rbind( DC_chrom_data, CC_chrom_data)
CC_EI_chrom_data <- rbind( EI_chrom_data, CC_chrom_data)

#### Gene Location Data ####

# to plot the links between the gene locations
# we need one dataframe for each genome
# MAKE SURE FILTERING IS THE SAME to ensure the rows match
# then calculate the midpoint of each gene and make a new df with chromosome and midpoint of the gene
Locations_Plot <- Locations_All_Plot %>% 
  filter(!is.na(CC_Chrom),
          CC_Chrom != 0 ) %>%
  rename( CC_start = start, CC_end = end ) %>% 
  mutate_at(vars(CM_start, CM_end, CC_start, CC_end, DC_start, DC_end), ~ str_replace_all(., ",", "" )) %>%
  mutate_at(vars(CM_start, CM_end, CC_start, CC_end, DC_start, DC_end), as.numeric ) %>%
  mutate( CM_mid = round( ( (CM_start) + (CM_end) )/2),
          CC_mid = round( ( (CC_start) + (CC_end) )/2),
          DC_mid = round( ( (DC_start) + (DC_end) )/2),
          EI_mid = round( ( (EI_start) + (EI_end) )/2),
          CM_Chrom = paste0("CM_", CM_Chrom),
          CC_Chrom = paste0("CC_", CC_Chrom),
          DC_Chrom = paste0("DC_", DC_Chrom),
          EI_Chrom = paste0("EI_", EI_Chrom))

# CM
CM_mid <- Locations_Plot %>% 
  select( CM_Chrom, CM_mid )

# CC
CC_mid <- Locations_Plot %>%
  select( CC_Chrom, CC_mid )

# DC
DC_mid <- Locations_Plot %>%
  select( DC_Chrom, DC_mid )

# EI
EI_mid <- Locations_Plot %>%
  select( EI_Chrom, EI_mid)

cols <- plasma(28)

#####################
#### Make Plots #####
#####################

## Makes the plots one at a time, not looped. Sorry.

# Save the plot to a PDF file
pdf(file.path(Plot_DIR, "Circos_Plot_CC_CM.pdf"), width = 10, height = 10) # Specify the file name and dimensions

###### CC CM
# prepare cell
circos.clear()
# the gap after puts a 1 degree gap after each chromosome apart from the 28th and 56th (last of each species), which have a bigger gap
circos.par(gap.after = c(rep(1, 27), 5, rep(1, 27), 5),
           "start.degree" = 87.5) # start degree puts first chromosomes at the top

# this intialises a blank plot, but the plot has dimensions according to the data provided
circos.genomicInitialize(CC_CM_chrom_data,
                         plotType = NULL # whether to plot outside of the chromosome with size indicators (Mb ticks)
)

circos.track(ylim = c(0, 1),
             bg.col = c(cols, rev(cols)),
             bg.border = NA, track.height = 0.05)

chromosome_to_color <- setNames(cols, mixedsort( unique(CC_mid$CC_Chrom)) )
color_vector <- chromosome_to_color[CC_mid$CC_Chrom] %>% as.vector()

# which genes are not syntenic between green and loggerhead
CC_CM_NonSynt_IDX <- which( gsub("CM_", "", CM_mid$CM_Chrom) != gsub("CC_", "", CC_mid$CC_Chrom))

# plot the normal lines
circos.genomicLink( CM_mid[ -CC_CM_NonSynt_IDX, ], CC_mid[ -CC_CM_NonSynt_IDX, ], col = color_vector[ -CC_CM_NonSynt_IDX]  )

# plot non-syntenic lines again, in dashed red, to emphasise 
circos.genomicLink( CM_mid[ CC_CM_NonSynt_IDX, ], CC_mid[ CC_CM_NonSynt_IDX, ], col = "red", lwd = 2, lty = "dashed"  )

text(1,0.8,"CM")
text(-1,0.8,"CC")

# Close the PDF device
dev.off()

##################

##### CC DC

# Save the plot to a PDF file
pdf(file.path(Plot_DIR, "Circos_Plot_CC_DC.pdf"), width = 10, height = 10) # Specify the file name and dimensions

# prepare cell
circos.clear()
# the gap after puts a 1 degree gap after each chromosome apart from the 28th and 56th (last of each species), which have a bigger gap
circos.par(gap.after = c(rep(1, 27), 5, rep(1, 27), 5),
           "start.degree" = 87.5) # start degree puts first chromosomes at the top

# this intialises a blank plot, but the plot has dimensions according to the data provided
circos.genomicInitialize(CC_DC_chrom_data,
                         plotType = NULL # whether to plot outside of the chromosome with size indicators (Mb ticks)
)

circos.track(ylim = c(0, 1),
             bg.col = c(cols, rev(cols)),
             bg.border = NA, track.height = 0.05)

chromosome_to_color <- setNames(cols, mixedsort( unique(CC_mid$CC_Chrom)) )
color_vector <- chromosome_to_color[CC_mid$CC_Chrom] %>% as.vector()

# which genes are not syntenic between green and loggerhead
CC_DC_NonSynt_IDX <- which( gsub("DC_", "", DC_mid$DC_Chrom) != gsub("CC_", "", CC_mid$CC_Chrom))

# plot the normal lines
circos.genomicLink( DC_mid[ -CC_DC_NonSynt_IDX, ], CC_mid[ -CC_DC_NonSynt_IDX, ], col = color_vector[ -CC_DC_NonSynt_IDX]  )

# plot non-syntenic lines again, in dashed red, to emphasise 
circos.genomicLink( DC_mid[ CC_DC_NonSynt_IDX, ], CC_mid[ CC_DC_NonSynt_IDX, ], col = "red", lwd = 2, lty = "dashed"  )

text(1,0.8,"DC")
text(-1,0.8,"CC")

# Close the PDF device
dev.off()

##### CC EI

# Save the plot to a PDF file
pdf(file.path(Plot_DIR, "Circos_Plot_CC_EI.pdf"), width = 10, height = 10) # Specify the file name and dimensions

# prepare cell
circos.clear()
# the gap after puts a 1 degree gap after each chromosome apart from the 28th and 56th (last of each species), which have a bigger gap
circos.par(gap.after = c(rep(1, 27), 5, rep(1, 27), 5),
           "start.degree" = 87.5) # start degree puts first chromosomes at the top

# this intialises a blank plot, but the plot has dimensions according to the data provided
circos.genomicInitialize(CC_EI_chrom_data,
                         plotType = NULL # whether to plot outside of the chromosome with size indicators (Mb ticks)
)


circos.track(ylim = c(0, 1),
             bg.col = c(cols, rev(cols)),
             bg.border = NA, track.height = 0.05)

chromosome_to_color <- setNames(cols, mixedsort( unique(CC_mid$CC_Chrom)) )
color_vector <- chromosome_to_color[CC_mid$CC_Chrom] %>% as.vector()

# which genes are not syntenic between green and loggerhead
CC_EI_NonSynt_IDX <- which( gsub("EI_", "", EI_mid$EI_Chrom) != gsub("CC_", "", CC_mid$CC_Chrom))

# plot the normal lines
circos.genomicLink( EI_mid[ -CC_EI_NonSynt_IDX, ], CC_mid[ -CC_EI_NonSynt_IDX, ], col = color_vector[ -CC_EI_NonSynt_IDX]  )

# plot non-syntenic lines again, in dashed red, to emphasise 
# note this is none for EI CC
#circos.genomicLink( EI_mid[ CC_EI_NonSynt_IDX, ], CC_mid[ CC_EI_NonSynt_IDX, ], col = "red", lwd = 2, lty = "dashed"  )

text(1,0.8,"EI")
text(-1,0.8,"CC")

# Close the PDF device
dev.off()
