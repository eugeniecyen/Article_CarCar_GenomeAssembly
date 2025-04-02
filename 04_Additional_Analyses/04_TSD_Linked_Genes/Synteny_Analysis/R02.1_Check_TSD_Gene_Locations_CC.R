### R02.2_Check_TSD_Gene_Locations_CC.R ###

########################################################################################################################################

# Created by: James Gilbert (https://github.com/evoJDG), 2024

### We want to compare in different species' genomes the locations of genes linked to sex determination
### In a previous step, we look up ~200 genes in the genomes of loggerhead (CC) and hawsbill (EI) turtles
### We already have data for green (CM) and leatherback (DC) turtles
### We looked up the genes by 1. BLASTing gene sequences against the genomes and 2) searching by name in the annotation file

### Here, we look through the CC data to verify 1) that the sequences we extracted are for the correct gene and
### 2) any instances of non-synteny are not due to the way we have found sequences.
# 1. Check genes by name
# 2. Merge with locations data from other species
# 3. Check non-syntenic
### At the end, we save to file the locations of SD genes for use in a later analysis.

#### Note, we do the same thing for the EI genes in a separate script. Is too messy to do both here.
#### This work is in response to reviewer comments.

########################################################################################################################################

###########################
##### Prep Environment ####
###########################

library(dplyr)

# Load files
## main project dir for Submission-2, the work we are doing in response to reviewers
DIR <- "/data/SBCS-EizaguirreLab/James/SD_Flex/Submission-2"

#### Locations of the CC SD genes, identified in previous step
SD_Locations_CC <- read.table(file.path(DIR, "01_Synteny/data/Locations/SDGene_Locations_CC.txt"), sep = "\t", header=T) %>%
  filter( !grepl("_contig_", chromosome)) %>%
  mutate(length = end - start + 1) # note start and end are inclusive, so we add 1

#### read in locations of genes in green and leatherback genomes
SD_locations_CM_DC <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/Bentley_TSDgene_info.txt"), sep = "\t", header = T) %>%
  # remove messy columns
  select(-starts_with("NUCL_SEQ"), -starts_with("PROT_SEQ"), -Reference) %>%
  # find duplicated genes in CM and DC
  group_by(GENE_ID) %>%
  mutate(duplicated_in_CM = ifelse(n_distinct(CM_Chrom) > 1, "Yes", "No"),
         duplicated_in_DC = ifelse(n_distinct(DC_Chrom) > 1, "Yes", "No")) %>%
  ungroup

################################################################
########## 1. Check Sequences Come From Correct Gene ###########
################################################################
## The format is a bit unusual
## The script reads in the data and throughout the script we check the data
## If we spot something that should be removed or confirmed as correct,
## we return to this point and add the info to an Action column.

# first look at genes whose query name (gene) doesn't matcht the name of the sequence extracted (gene_extracted)
SD_Locations_CC %>% filter(genes_match == "different") %>%
  mutate(gene = as.factor(gene)) %>%
  select(gene, sequence, gene_extracted, genes_match)

######### Define Actions #########
######### this is an important bit

### Adding remove to the action var means they are removed in a later step
### i.e. this is a way of removing incorrect sequences or genes

## If gene names don't match, we look up the details in the annotation
## And manually blast to see whether the sequence matches the target gene
## where the sequence doesn't correspond to the query gene, we flag it for removal

## we also check genes where the chromosomes are different between CC and CM and/or DC
## this might be because the sequence extracted may not correspond to the query gene
## e.g. it might be a different gene in the same family located on a different chromosome, or it might be a different gene entirely
## We check non-syntenic cases to avoid false positive identification of non-synteny
## For non-synteny cases, we check the extracted sequence in annotation and by BLASTing,
## to confirm whether the sequence corresponds to the query gene

SD_Locations_CC <- SD_Locations_CC %>% mutate( action = NA,
  action = case_when( 
  
  #### gene IDs don't match ####
  gene == "CYP19A" & sequence == "CarCarScaffG00000020878.1" ~ "Correct gene. No action",
  gene == "DAX1" & gene_extracted == "NR0B1" ~ "Aliases. No action",
  gene == "DAZL" & sequence == "CarCarScaffG00000005519.1" ~ "Wrong gene. Remove",
  gene == "EPOR" & gene_extracted == "RGL3" ~ "Different genes. Remove",
  gene == "FAM132A" & gene_extracted == "UBE2J2" ~ "Different genes. Remove",
  gene == "ELMO1" & sequence == "CarCarScaffG00000007577.1" ~ "Blast returns diff genes. Remove",
  gene == "HEMGN" & sequence == "CarCarScaffG00000013150.1" ~ "Blast confirms correct gene. No action",
  gene == "HSPB6" & gene_extracted == "CRYAB" ~ "Different genes. Remove",
  gene == "MMP23B" & gene_extracted == "MMP23" ~ "Aliases. No action",
  gene == "MYL1" & gene_extracted == "RTASE" ~ "Different genes. Remove",
  gene == "UPK3A" & gene_extracted == "FAM118A" ~ "Different genes. Remove",
  
  ####  chromosomes don't match with CM/DC   #### 
  gene == "ACSS1" & sequence == "CarCarScaffG00000014799.1" ~ "Seq annotated as Similar to ACSS1 but BLAST confirms ACSS2. Remove.",
  gene == "EP300" & chromosome == "SLK063_ragtag_chr1" ~ "Gene mapped to two locations in CM and DC but only once in CC. No action required.",
  gene == "SLC23A1" & sequence == "CarCarScaffG00000003657.1" ~  "Seq annotated as Similar to SLC23A1 Blast returns C mydas SLC23A1. Looks like correct gene but might be incorrectly placed - Remove",
  # SLC23A1 is removed following logic outlined in the methods - where there are hits on 2 chroms in CC but only 1 in CM/DC, we remove the non-syntenic hit.
  # Conservatively assuming annotation error.
  gene == "HSPA8"  ~  "Both versions in CM/DC produce hits in CC. Looks like genuine duplication. No action",
  gene == "GSK3B" & sequence == "CarCarScaffG00000029466.1" ~ "Blast returns  GSK3A. Remove",
  gene == "MAP3K3" & sequence %in% c("CarCarScaffG00000004917.1",
                                     "CarCarScaffG00000004918.1",
                                     "CarCarScaffG00000004919.1",
                                     "CarCarScaffG00000030860.1") ~ "Gene mapped to two locations in all species. Looks like genuine duplication. No action required.",
  gene == "CIRBP" & sequence %in% c("CarCarScaffG00000031169.1",
                                    "CarCarScaffG00000031170.1",
                                    "CarCarScaffG00000031171.1") ~  "Seq annotated as Similar to CIRBP Blast returns C mydas CIRBP. Correct gene, no action",
  gene == "IFIT5" & sequence == "CarCarScaffG00000016829.1" ~  "Seq annotated as Similar to IFIT5 Blast returns C mydas IFIT5. No action"
  ))

# Remove sequences that do not match with query gene
SD_Carcar <- SD_Locations_CC %>% filter( !grepl("Remove", action))

# Now we have removed incorrect sequences, we don't need multiple sequences (i.e. isoforms) per gene
# select the longest seq per gene
SD_Carcar <- SD_Carcar %>% 
  group_by(gene, chromosome) %>%
  slice_max(order_by = length, n = 1, with_ties = T) %>% # Select longest seqs
  arrange(sequence) %>% # Ensure seq_id ordering for seqs of same length
  slice_head(n = 1) %>% # Take the first if lengths are equal
  ungroup()

## check whether gene duplicated (found on more than one chromosome - don't count isoforms as duplications)
SD_Carcar <- SD_Carcar %>%
  group_by(gene) %>%
  mutate(duplicated_in_CC = ifelse(n_distinct(chromosome) > 1, "Yes", "No") ) %>%
  ungroup

# Check genes whose names don't match query
SD_Carcar %>% filter(genes_match == "different") %>%
  select(gene, sequence, gene_extracted, genes_match, action)

#########################################################
########## 2. Map chromosomes between species ###########
#########################################################

# To know whether genes are syntenic, we need to know which chromosomes are equivalent in the different genomes
# The chromosomes have different names in the SD gene metadata (e.g. SUPER_1) and synteny analysis (NC_12345)

# To get from Super_1 to the corresponding loggerhead chromosome:
# 1. Map SUPER_1 (Bentley metadata) to Chromosome sequence ID (chromosome mapping file, NCBI)
# 2. Map sequence ID to corresponding loggerhead chromosome (synteny analysis file)
# 3. Add to SD seq metadata

##### For green turtle #####
# 1. chromosome mapping file
Chr_map_CM <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/CM_Chrom_Map.tsv"), sep = "\t", header = T) %>%
  select(Chromosome.name, RefSeq.seq.accession) %>%
  filter(grepl("NC_", RefSeq.seq.accession))

# 2. CC to CM chromosome mapping from synteny analysis
Chr_map_CM_CC <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/CheMyd_CarCar_assoc.tsv"), sep = "\t", header = T) %>%
  filter(grepl("NC_", Query)) %>% # keep only full chromosomes that have NC prefix
  select(Query, Target) %>% # remove unnecessary columns
  left_join(x = ., 
            y = Chr_map_CM,
            join_by(Query == RefSeq.seq.accession ) ) %>% # match correct columns
  select(Query, Target, Chromosome.name)


##### For leatherback turtle #####
# 1. chromosome mapping file
Chr_map_DC <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/DC_Chrom_Map.tsv"), sep = "\t", header = T) %>%
  select(Chromosome.name, GenBank.seq.accession) %>% # NOTE the different chromosome seq identifiers - CM uses RefSeq, DC uses GenBank
  filter(grepl("CM", GenBank.seq.accession))

# 2. CC to CM chromosome mapping from synteny analysis
Chr_map_DC_CC <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/DerCor_CarCar_assoc.tsv"), sep = "\t", header = T) %>%
  filter(grepl("CM", Query)) %>% # keep only full chromosomes that have CM prefix
  select(Query, Target) %>%
  left_join(x = ., 
            y = Chr_map_DC,
            join_by(Query == GenBank.seq.accession ) ) %>%
  select(Query, Target, Chromosome.name)

##### Combine Locations for All Species #####
## the important thing is that for each chromosome in the SD locations for CM and DC,
## we have the equivalent chromosome in CC. This is what we compare for synteny
# 3. Add CC equivalent in SD seq metadata
SD_locations_CM_DC = SD_locations_CM_DC %>%
  mutate(CM_Chrom_cleaned = gsub("SUPER_", "", CM_Chrom),
         DC_Chrom_cleaned = gsub("SUPER_", "", DC_Chrom))  %>%
  ## add chromosome equivalent for CM
  left_join(
    Chr_map_CM_CC %>% # add map CM CC
      select(Chromosome.name, Target), # only need chrm number and CC Target
    by = c("CM_Chrom_cleaned" = "Chromosome.name") # Match cleaned CM_Chrom to Chromosome.name
         ) %>%
  rename( CM_CC_Equivalent = Target) %>%
  as.data.frame

SD_locations_CM_DC = SD_locations_CM_DC %>%
  mutate(DC_Chrom_cleaned = as.integer(DC_Chrom_cleaned)) %>%
  ## add chromosome equivalent for DC
  left_join(
    Chr_map_DC_CC %>% # add map DC CC
      select(Chromosome.name, Target), # only need chrm number and CC Target
    by = c("DC_Chrom_cleaned" = "Chromosome.name") # Match cleaned CM_Chrom to Chromosome.name
  ) %>%
  rename( DC_CC_Equivalent = Target) %>%
  as.data.frame %>%
  select(-CM_Chrom_cleaned, -DC_Chrom_cleaned)

### add CC locations to main data
SD_Locations_All <- merge(x = SD_Carcar, by.x = "gene",
                       y = SD_locations_CM_DC %>% select(c(GENE_ID, starts_with("DC_"), starts_with("CM_"), duplicated_in_CM, duplicated_in_DC) ), by.y = "GENE_ID",
                       all.x = T, all.y = T) %>%
  select(gene, sequence, chromosome, start, end, starts_with("CM_"), starts_with("DC_"), BlastOrName, action, duplicated_in_CC, duplicated_in_CM, duplicated_in_DC) %>%
  mutate( CC_Chrom = as.numeric(gsub("SLK063_ragtag_chr", "", chromosome)),
          CM_Chrom = as.numeric(gsub("SLK063_ragtag_chr", "", CM_CC_Equivalent)),
          DC_Chrom = as.numeric(gsub("SLK063_ragtag_chr", "", DC_CC_Equivalent)),
          found_in_CC = if_else(is.na(sequence), "No", "Yes"),
          match_CC_CM = if_else(chromosome == CM_CC_Equivalent, "Yes", "No"),
          match_CC_DC = if_else(chromosome == DC_CC_Equivalent, "Yes", "No")) %>% unique

######################################################
########## 3. Check non-syntenic sequences ###########
######################################################

#### Make sure sequences that are not syntenic represent correct genes ####
### don't want false positive caused by a non-syntenic sequence assigned to the wrong gene

### 1. Look at the genes that are non syntenic
### 2. Look up info in annotation file
### 3. And blast sequences when not sure
### If seq does not correspond to correct gene, add to remove section at the start of the script
### If it is the correct gene and looks like genuine non-synteny or duplication, also add a note above
SD_Locations_All %>% filter(match_CC_CM == "No")
SD_Locations_All %>% filter(match_CC_DC == "No")

### HSPA8 and MAP3K3 are duplicated in all species.
### each sequence i.e. CC_HSPA8_1 maps to both duplicates in the other species
### so you get CC_HSPA8_1 = CM_HSPA8_1 and CM_HSPA8_2
### and CC_HSPA8_2 = CM_HSPA8_1 and CM_HSPA8_2
### need to remove the non-syntenics:
REMOVE <- which( SD_Locations_All$gene %in% c("HSPA8","MAP3K3") & SD_Locations_All$match_CC_CM == "No" )

# check these should be removed:
SD_Locations_All[ REMOVE, ]

# remove em
SD_Locations_All <- SD_Locations_All [ -REMOVE, ]

# does it look sensible
## overview
SD_Locations_All %>% select( gene, ends_with("Chrom")) %>% arrange(CC_Chrom)
## duplications in CC
SD_Locations_All %>% filter(duplicated_in_CC == "Yes") %>% select(gene, sequence, chromosome, start, end, duplicated_in_CC, match_CC_CM, match_CC_DC)
## non-syntenic cases
SD_Locations_All %>% filter(chromosome != CM_CC_Equivalent)

SD_Locations_All %>% filter(chromosome != DC_CC_Equivalent)

# save output
write.table(SD_Locations_All, file.path(DIR, "01_Synteny/data/Locations/SDGene_Locations.out.txt"), sep = "\t", row.names = F)

#####################################################
############# Locations for Methylation #############
#####################################################

# We will later look at the methylation of SD genes in CC.
# only need the locations from the CC genome
# and for genes that are duplicated in CM/DC but not CC, remove one of the rows so the same gene (and sequence) is not included twice in methylation
SDLocations_forMethylation = SD_Locations_All %>% filter( found_in_CC == "Yes" ) %>%
  select(gene, sequence, chromosome, start, end) %>% unique


# Check the number of rows in the data frame = 201
if (nrow(SDLocations_forMethylation) != 201) {
  stop("Error: The number of rows in SDLocations_forMethylation is not 201.")
}

## should only by HSPA8 or MAP3K3 - genuine duplicates
SDLocations_forMethylation %>% group_by(gene) %>% summarise(n = n()) %>% filter(n > 1)


### save to file
write.table(SDLocations_forMethylation,
            file.path(DIR, "01_Synteny/data/Locations/SDGene_Locations.forMethylationAnalysis.txt"), sep = "\t", row.names = F)



