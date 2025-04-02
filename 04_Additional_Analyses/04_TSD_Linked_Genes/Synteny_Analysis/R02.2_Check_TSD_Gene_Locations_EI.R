### R02.2_Check_TSD_Gene_Locations_EI.R ###

########################################################################################################################################

# Created by: James Gilbert (https://github.com/evoJDG), 2024

### Following S01 and S02 scripts, here we check the sequences and locations of genes in the hawksbill genome
### Full details are given in S02, where we do the same for the loggerhead genome
### In short, we extracted genes that matched name or sequence (BLAST)
### with a list of genes with documented links to sex determination
### Here, we check the sequences correspond to the intended gene
### Map the locations of the genes between CC and EI genomes
### Check any non-sytenic cases are not caused by the process of finding them

########################################################################################################################################

###########################
##### Prep Environment ####
###########################

library(dplyr)

# Load data
DIR <- "/data/SBCS-EizaguirreLab/James/SD_Flex/Submission-2"

### locations of genes in EI genome (identified in S01)
SD_Locations_EI <- read.table(file.path(DIR, "01_Synteny/data/Locations/SDGene_Locations_EI.txt"), sep = "\t", header=T) %>%
  mutate(length = end - start + 1) %>% # note start and end are inclusive, so we add 1
  rename(EI_start = start,
         EI_end = end) # to avoid clashing with var names in CC locations

### Locations of SD genes in all species, compiled in S02
SD_Locations_All = read.csv(file.path(DIR, "01_Synteny/data/Locations/SDGene_Locations.out.txt"), sep = "\t" )

################################################################
########## 1. Check Sequences Come From Correct Gene ###########
################################################################
## The format is a bit unusual
## The script reads in the data and throughout the script we check the data
## If we spot something that should be removed or confirmed as correct,
## we return to this point and add the info to an Action column.

# first look at genes whose query name (gene) doesn't matcht the name of the sequence extracted (gene_extracted)
SD_Locations_EI %>% filter(genes_match == "different") %>%
  mutate(gene = as.factor(gene)) %>%
  select(gene, sequence, gene_extracted, genes_match)

### note gene annotations hawksbill are not in GFF but in $EI_DIR/NSCA_data/Eretmochelys_imbricata.gene.function.xls
SD_Locations_EI <- SD_Locations_EI %>% mutate( action = NA,
                                               action = case_when( 
                                                 # gene IDs don't match
                                                 
                                                 ### sequences that produce a lot of hits - remove means no wrong hits in the data
                                                 ### but not sure why they produce so many hits. Verified it happens both in our analysis and on BLAST website:
                                                 ### They do hit many genes
                                                 gene != "ACAA2" & sequence == "Eim15045.1" ~
                                                   "This seq is best hit for 2 genes that are not ACAA2 - manual blast also produces lots of hits. Remove.",
                                                # sequence == "Eim13399.1" ~ "Sequence produces many hits in EI genome and on BLAST website - Remove.",
                                                 gene != "KATNB1" & sequence == "Eim04266.1" ~ "Sequence produces many hits in EI genome and on BLAST website - Remove.",
                                                 
                                                 ## Genes that don't match
                                                 gene == "ACSS1" & sequence == "Eim13399.1" ~ "BLAST returns uncharacterised, too. Can't verify so remove.",
                                                 gene == "ANPEP" & sequence == "Eim04266.1" ~ "Wrong gene. Remove.",
                                                 gene == "AXIN2" & sequence == "Eim05335.1" ~ "Same gene. No action.",
                                                 gene == "CLTRN" & sequence == "Eim00958.1" ~ "Same gene. No action",
                                                 gene == "COL11A1" & sequence == "Eim03341.1" ~ "Wrong gene. Remove",
                                                 gene == "CRABP2" & sequence == "Eim11176.2" ~ "Manual BLAST confirms Wrong gene. Remove.",
                                                 gene == "CYP19A" & sequence == "Eim03344.1" ~ "Same gene. No action",
                                                 gene == "DAX1" & sequence == "Eim01010.1" ~ "Aliases, no action",
                                                 gene == "DTNA" & sequence == "Eim17152.1" ~ "Manual BLAST suggests Wrong gene. Remove.",
                                                 gene == "FANK1" & sequence == "Eim08133.1" ~ "Manual BLAST: can't confirm FANK1. Remove.",
                                                 gene == "FOXL2" & sequence == "Eim19601.1" ~ "Same gene. No action",
                                                 gene == "GGT1" & sequence == "Eim06320.1" ~ "Wrong gene. Remove.",
                                                 gene == "MYH11" & sequence == "Eim00976.1" ~ "Wrong gene. Remove.",
                                                 gene == "PDGFA" & sequence == "Eim15870.1" ~ "Low quality protein. Remove.",
                                                 gene == "TMEM35A" & sequence == "Eim19858.1" ~ "Same gene. No action",
                                                 gene == "TNNI2" & sequence == "Eim16329.1" ~ "Manual BLAST confirms Wrong gene. Remove.",
                                                 gene == "TSPAN9" & sequence == "Eim03341.1" ~ "Wrong gene. Remove.",
                                                 gene == "UPK3A" & sequence == "Eim01723.1" ~ "Manual BLAST suggests correct gene. No action",
                                                 
                                                ## Chromosomes that don't match
                                                ## note some of these are annotated differently in the EI annotation to what is returned by BLAST
                                                gene == "A2M" & sequence == "rna-Eim12149.1" ~ "manual BLAST confirms A2M-like, not A2m. Different gene, Remove",
                                                gene == "ACSS1" & sequence == "rna-Eim16820.1" ~ "ACSS2. Different gene, Remove",
                                                gene == "CADPS" & sequence == "rna-Eim02391.1" ~ "CADPS2. Different gene, Remove",
                                                gene == "COL11A1" & sequence == "rna-Eim09837.1" ~ "Manual BLAST confirms COL5A3. Different gene, Remove",
                                                gene == "COL4A3" & sequence == "rna-Eim15253.1" ~ "Manual BLAST confirms CERT1. Different gene, Remove",
                                                gene == "COL4A3" & sequence == "rna-Eim15253.2" ~ "Manual BLAST confirms CERT1. Different gene, Remove",
                                                gene == "COL9A1" & sequence == "rna-Eim08407.1" ~ "Manual BLAST confirms different gene, Remove",
                                                gene == "FOXL2" & sequence == "rna-Eim07305.1" ~ "Not FOXL2. FOXL2-like. Remove",
                                                gene == "FOXL2" & sequence == "rna-Eim11647.1" ~ "Not FOXL2. FOXL2-like. Remove",
                                                gene == "GSK3B" & grepl("^rna-Eim11022.*", sequence) ~ "STARD10, not STAR. Remove",                                                gene == "SLC23A1" & sequence == "rna-Eim02148.1" ~ "One of two locations for SLC23A1 in EI. Following logic set out in methods, Remove the non-syntenic location",
                                                gene == "STAR" & grepl("^rna-Eim00292.*", sequence) ~ "STARD10, not STAR. Remove",
                                                gene == "STAR" & grepl("^rna-Eim19910.*", sequence) ~ "STARD10, not STAR. Remove",
                                                gene == "STAR" & sequence == "rna-Eim12289.1" ~ "STARD3, not STAR. Remove",
                                                gene == "STAR" & sequence == "rna-Eim15362.1" ~ "STARD4, not STAR. Remove",
                                                gene == "STAR" & sequence == "rna-Eim16658.1" ~ "STARD9, not STAR. Remove",
                                                gene == "STAR" & sequence == "rna-Eim17969.1" ~ "STARD1, not STAR. Remove",
                                                gene == "STAR" & sequence == "rna-Eim18701.1" ~ "STARD10, not STAR. Remove",
                                                gene == "STAR" & sequence == "rna-Eim19924.1" ~ "STARD8, not STAR. Remove",
                                                gene == "STAR" & sequence == "rna-Eim19925.1" ~ "STARD8, not STAR. Remove",
                                                gene == "STAR" & sequence == "rna-Eim00559.1" ~ "STARD13, not STAR. Remove",
                                                
                                               ))

# write output after removing incorrect genes
SD_Ereimb <- SD_Locations_EI %>% filter( !grepl("Remove", action),
                                         !grepl("unanchor", chromosome)) %>% # don't want genes on unanchored contigs 
  unique() 

## check whether gene duplicated (found on more than one chromosome - don't count isoforms as duplications)
SD_Ereimb <- SD_Ereimb %>%
  group_by(gene, chromosome) %>%
  slice_max(order_by = length, n = 1, with_ties = T) %>% # Select longest seqs
  arrange(sequence) %>% # Ensure seq_id ordering for seqs of same length
  slice_head(n = 1) %>% # Take the first if lengths are equal
  ungroup()

## check whether gene duplicated (found on more than one chromosome - don't count isoforms as duplications)
SD_Ereimb <- SD_Ereimb %>%
  group_by(gene) %>%
  mutate(duplicated_in_EI = ifelse(n_distinct(chromosome) > 1, "Yes", "No") ) %>%
  ungroup

# Check genes whose names don't match query
SD_Ereimb %>% filter(genes_match == "different") %>%
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

# To know whether genes are syntenic, we need to know which chromosomes are equivalent in the different genomes

# 1. chromosome mapping file
Chr_map_EI <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/EI_Chrom_Map.tsv"), sep = "\t", header = T) %>%
  select(Chromosome.name, Sequence.name, GenBank.seq.accession) %>%
  filter(grepl("CM",  GenBank.seq.accession))

# CC to EI chromosome mapping from synteny analysis
Chr_map_EI_CC <- read.table( file.path(DIR, "Genome_SD_Flex/00_Metadata/EreImb_CarCar_assoc.tsv"), sep = "\t", header = T) %>%
  filter(grepl("CM", Query)) %>% # keep only full chromosomes that have CM prefix
  select(Query, Target) %>%
  merge(x = .,
            by.x = "Query",
            y = Chr_map_EI, by.y = "GenBank.seq.accession")


##### Combine Locations for CC and EI #####
## the important thing is that for each chromosome in the SD locations for EI,
## we have the equivalent chromosome in CC. This is what we compare for synteny

# 3. Add CC equivalent to EI seq data
SD_Ereimb = SD_Ereimb %>%
  left_join(
    Chr_map_EI_CC %>% select(Target, Sequence.name),
    # CC chrom is called chromosome in SD_Locations_All and Target in the EI CC map
    by = c("chromosome" = "Sequence.name")
  ) %>%
  rename( EI_CC_Equivalent = Target) %>%
  as.data.frame


### add EI locations to main data
SD_Locations_All_EI = merge(x = SD_Locations_All %>% select(gene, chromosome) %>% unique,
                         y = SD_Ereimb,
                         by = "gene") %>% 
  rename(chromosome_CC = chromosome.x,
         chromosome_EI = chromosome.y) %>%
  mutate( CC_Chrom = as.numeric(gsub("SLK063_ragtag_chr", "", chromosome_CC)),
          EI_Chrom = as.numeric(gsub("SLK063_ragtag_chr", "", EI_CC_Equivalent)),
          found_in_EI = if_else(is.na(sequence), "No", "Yes"),
          match_CC_EI = if_else(chromosome_CC == EI_CC_Equivalent, "Yes", "No")) %>%
  unique %>%
  arrange(gene)

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
SD_Locations_All_EI %>% filter(match_CC_EI == "No")

### for CC EI it is simpler than CM and DC. The only non-matches are HSPA8 and MAP3K3.
### Both genes are duplicated and the non-match is between paralogues on different chromosomes
### Checked annotation and both genes are duplicated and present on same chrms in both species:

## HSPA8 is on 21 and 22 in CC. EI blast returns Eim11508 and Eim17031 and Eim10261.
## 11508 is HSPA8 (synonym in annotation - not picked up in BLAST results) on Superscaffold 24 (SLK063_ragtag_chr22).
## 17031 is HSPA2 i.e. different gene
## 10261 is HSPA8 on Superscaffold 22 (SLK063_ragtag_chr21)

## MAP3K3 on 2 and 24 in CC. EI blast returns one seq: Eim12174.1 (Superscaffold2 = SLK063_ragtag_chr2).
## Also in annotation Eim12174 (Superscaffold 27 = SLK063_ragtag_chr24)

### Given I manually checked these are the same in both species, we will move the non-matching sequences here.
REMOVE <- which( SD_Locations_All_EI$gene %in% c("HSPA8","MAP3K3") & SD_Locations_All_EI$match_CC_EI == "No" )

# check these should be removed:
SD_Locations_All_EI[ REMOVE, ]

# remove em
SD_Locations_All_EI <- SD_Locations_All_EI [ -REMOVE, ]

# does it look sensible
## overview
SD_Locations_All_EI %>% select( gene, ends_with("Chrom")) %>% arrange(CC_Chrom)
## number found in SD_Locations_All_EI
SD_Locations_All_EI %>% select( gene, ends_with("Chrom")) %>% nrow()
## duplications in CC
SD_Locations_All_EI %>% filter(duplicated_in_EI == "Yes") %>% select(gene, sequence, chromosome_EI)
## non-syntenic cases
SD_Locations_All_EI %>% filter(match_CC_EI == "No")

### make a df for the plot
### need for each gene: gene, chromosome_CC, chromosome_EI, EI_start, EI_end
### plus Locations_All
SD_Locations_All_Plot = merge( SD_Locations_All,
                               SD_Locations_All_EI,
                               by = c("gene", "CC_Chrom"),
                               all.x = T) %>%

  select(gene, ends_with("_Chrom"), ends_with("start"), ends_with("end"))

# save output
write.table(SD_Locations_All_Plot, file.path(DIR, "01_Synteny/data/Locations/SDGene_Locations_All_Plot.txt"), sep = "\t", row.names = F)



