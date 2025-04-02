### R03_Compare_ONT_WGBS_Methylomes.R ###

########################################################################################################################################

# Created by Alice, 2024 (https://github.com/alicebalard); edited by Charley, 2024

# Perform comparisons of methyation per CpG site and per gene by feature type between the ONT reference methylome
# and average WGBS methylome from 10 additional nesting females (including separate comparisons per WGBS individual)

########################################################################################################################################

###############################
###### Prep environment #######
###############################

library(data.table)
library(ggplot2)
library(viridis)

setwd("/data/SBCS-EizaguirreLab/")

########################################################################################################################################

###########################
###### Prep objects #######
###########################

rerunCalc = FALSE
if (rerunCalc){
  
  ########################
  ## Read in ONT methylome
  
  ## 1. Reading the methylation calls from sorted Guppy alignments (should be bam/sam like Bismark)
  ## Use fread, faster than read.csv (NB: tab delimited in bash)
  load("/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/myTurtleBedSum.RData")
  
  ## modkit pileup: 0-based start position https://github.com/nanoporetech/modkit?tab=readme-ov-file#description-of-bedmethyl-output
  ## need to add +1 to match Bismark 1-based convention
  myTurtleBedSum$pos = paste0(myTurtleBedSum$chrom, "_", myTurtleBedSum$start + 1)
  
  range(myTurtleBedSum$coverage)
  
  ######################################
  ## Ensure min. coverage of 8 or more ##
  myTurtleBedSum <- myTurtleBedSum[myTurtleBedSum$coverage >= 8,]
  
  ## Calculate percentage methylation to compare with WGBS
  myTurtleBedSum$percMeth = myTurtleBedSum$fracMeth * 100
  
  ## Load feature information
  DIR_ONT <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/01_Working/02_ONT"
  ONT_GRangeOBJ <- readRDS(file = file.path(DIR_ONT, "GRanges_myTurtleBedSum_GeneInfoPerSite.RDS"))
  ONT_features = data.frame(pos = paste0(ONT_GRangeOBJ@seqnames, "_", ONT_GRangeOBJ@ranges@start+1), ## add +1 
                            featureType = ONT_GRangeOBJ$featureType, geneInfo = ONT_GRangeOBJ$geneInfo)
  
  ONT_full = merge(data.frame(pos = myTurtleBedSum$pos, percMeth = myTurtleBedSum$percMeth),
                   ONT_features)
  
  saveRDS(ONT_full, 
          file = "/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/ONT_full.RDS")
  
  ############################
  ## Read in WGBS methylome data
  DIR_UniteCov <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/01_Working/03_WGBS/02_UniteCov_Objects"
  uniteCov10ind <- readRDS(file.path(DIR_UniteCov, "GenomePaper_WholeGenome_uniteCov75pc_Reloc3_Mothers.RDS")) 
  
  range(uniteCov10ind$coverage1, na.rm = T)
  
  ######################################
  ## Ensure coverage of 8 or more to match ONT ##
  df = getData(uniteCov10ind)
  
  # Create a list of column prefixes
  prefixes <- c("coverage", "numCs", "numTs")
  
  # Loop through each individual
  for (i in 1:10) {
    # Create column names
    cols <- paste0(prefixes, i)
    # Check if coverage column exists
    if (cols[1] %in% colnames(df)) {
      # Create a logical vector for coverage < 8
      low_coverage <- df[[cols[1]]] < 8 & !is.na(df[[cols[1]]])      
      # Replace values with NA where coverage is < 8
      df[low_coverage, cols] <- NA
    }
  }
  
  ## Replace in the methylobject
  for (i in 1:3){
    for (j in 1:10)
      uniteCov10ind[[paste0(prefixes[i], j)]]=df[[paste0(prefixes[i], j)]]
  }
  
  ## Calculate percentage of methylation
  WGBS_percMeth <- methylKit::percMethylation(uniteCov10ind)
  
  ## Add the positions names
  rownames(WGBS_percMeth) <- paste0(uniteCov10ind$chr, "_", uniteCov10ind$start)
  
  ## Filtrate for NA in 3 individuals maximum
  WGBS_percMeth = WGBS_percMeth[rowSums(!is.na(WGBS_percMeth)) >= 7 , ]
  
  ## Load feature information
  DIR_GRanges <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/01_Working/03_WGBS/03_GRanges_Objects"
  WGBS_GRangeOBJ <- readRDS(file = file.path(DIR_GRanges, "GRanges_75pc_Reloc3_Mothers_GeneInfoPerSite.RDS"))
  
  WGBS_features = data.frame(pos = paste0(WGBS_GRangeOBJ@seqnames, "_", WGBS_GRangeOBJ@ranges@start), 
                             featureType = WGBS_GRangeOBJ$featureType, geneInfo = WGBS_GRangeOBJ$geneInfo)
  
  WGBS_percMeth=data.frame(WGBS_percMeth)
  WGBS_percMeth$pos = rownames(WGBS_percMeth)
  
  WGBS_full = merge(WGBS_percMeth, WGBS_features)
  
  saveRDS(WGBS_full, 
          file = "/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/WGBS_full.RDS")
  
  ####################
  ## Calculate overlap
  getOverlap <- function(WGBSpos, ONTpos){
    totalSeq <- length(unique(c(WGBSpos, ONTpos))) 
    overlap <- length(intersect(WGBSpos, ONTpos)) 
    return(paste0("TotalSeq=", totalSeq, ", overlap=", overlap, ", percOverlap=", 
                  round(overlap/totalSeq*100, 3)))
  }
  
  ## ONT vs all WGBS positions sequenced in 10 individuals (cov >=8)
  # (BEFORE filtering for 3 NAs or less)
  posSequencedONT=myTurtleBedSum$pos
  
  posSequencedWGBS = paste0(uniteCov10ind$chr, "_", uniteCov10ind$start)[
    rowSums(is.na(getData(uniteCov10ind)[paste0("coverage", 1:10)])) != 10]
  
  getOverlap(posSequencedONT, posSequencedWGBS)
  # [1] "TotalSeq=26368272, overlap=23806545, percOverlap=90.285"
  
  ## ONT vs WGBS positions sequenced in each individual (cov >=8)
  for (i in 1:10){
  posWGBS=paste0(uniteCov10ind$chr, "_", uniteCov10ind$start)[
     !is.na(getData(uniteCov10ind)[paste0("coverage", i)])]
  print(i)
  print(getOverlap(posSequencedONT, posWGBS))
  }
     
  ## Individuals:
  uniteCov10ind@sample.ids
    
  rm(WGBS_percMeth, WGBS_features, WGBS_GRangeOBJ, uniteCov10ind, 
     myTurtleBedSum,ONT_GRangeOBJ,ONT_features)
}

########################################################################################################################################

###########################
###### Per CpG site #######
###########################

######## ONT vs population average WGBS methylome ########

### Prep data ###

# # DF for 1 individual ONT, with feature annotation:
# ONT_full = readRDS("/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/ONT_full.RDS")
# # DF for 10 individuals WGBS, positions sequenced in at least 7 (i.e. 75% as integer), with feature annotation:
# WGBS_full = readRDS("/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/WGBS_full.RDS")

# Calculate mean across all 10 individuals per CpG site
WGBS_full_average10 = data.frame(pos = WGBS_full$pos, 
                                 featureType = WGBS_full$featureType,
                                 geneInfo = WGBS_full$geneInfo,
                                 percMeth10ave = rowMeans(WGBS_full[,grep("SLL", names(WGBS_full))], na.rm = T))

combined_ONT_WGBS10ave = merge(ONT_full, WGBS_full_average10)

### Linear model ###

lm <- lm(percMeth ~ percMeth10ave*featureType, data=combined_ONT_WGBS10ave)
anova(lm)

### Pearson's test ###

# -> Interaction by feature type -> Pearson's correlation separately per feature type as a post-hoc test
# Subset
exons <- combined_ONT_WGBS10ave[combined_ONT_WGBS10ave$featureType == "exons",]
introns <- combined_ONT_WGBS10ave[combined_ONT_WGBS10ave$featureType == "introns",]
promoters <- combined_ONT_WGBS10ave[combined_ONT_WGBS10ave$featureType == "promoters",]
intergenic <- combined_ONT_WGBS10ave[combined_ONT_WGBS10ave$featureType == "intergenic",]

# exons
cor_exon <- cor.test(exons$percMeth, exons$percMeth10ave)
data.frame(t=cor_exon$statistic, df=cor_exon$parameter, 
           p.value=cor_exon$p.value, R2=round(cor_exon$estimate, 3))

# introns
cor_intron <- cor.test(introns$percMeth, introns$percMeth10ave)
data.frame(t=cor_intron$statistic, df=cor_intron$parameter, 
           p.value=cor_intron$p.value, R2=round(cor_intron$estimate, 3))

# promoters
cor_prom <- cor.test(promoters$percMeth, promoters$percMeth10ave)
data.frame(t=cor_prom$statistic, df=cor_prom$parameter, 
           p.value=cor_prom$p.value, R2=round(cor_prom$estimate, 3))

# intergenic
cor_int <- cor.test(intergenic$percMeth, intergenic$percMeth10ave)
data.frame(t=cor_int$statistic, df=cor_int$parameter, 
           p.value=cor_int$p.value, R2=round(cor_int$estimate, 3))

### Plot ###
viridis_custom_palette <- c("#0d0887", "#5b01a5", "#a82296", "#da5b69", "#fa9e3b")
OUTDIR_Figures <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/02_Output/Figures"

# Randomly subset 100,000 sites and plotting 
set.seed(1234)
scatterPlot_byCpG <- ggplot(combined_ONT_WGBS10ave[sample(1:nrow(combined_ONT_WGBS10ave), 100000),], 
                            aes(x=percMeth, y=percMeth10ave, color=featureType) ) +
  geom_point(alpha=0.5, size=0.1) +
  stat_smooth(method = "lm", color="black", alpha=0.8, linewidth = 0.8) +
  labs(x="ONT: methylation (%)", y="WGBS: mean methylation for 10 individuals (%)", color="Feature type") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  theme_classic() +
  scale_color_manual(values = c(viridis_custom_palette[5], viridis_custom_palette[1], viridis_custom_palette[2], viridis_custom_palette[3])) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16), legend.position="none") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  facet_wrap(~featureType, 
             labeller = labeller(featureType = c("exons" = "Exon", "intergenic" = "Intergenic (<10kb from TSS)", "introns" = "Intron", "promoters" = "Promoter"))) +
  theme(strip.text.x = element_text(size = 14))

ggsave(file=file.path(OUTDIR_Figures ,"By_Site/ScatterPlot_ONT_vs_WGBS_MeanMeth_AnnotCpGs_ByFeatureType_100ksites.pdf"), 
       plot=scatterPlot_byCpG , 
       width=7, height=7, units="in")

######## Individually per WGBS methylome ########

combined_ONT_WGBS = merge(ONT_full, WGBS_full)
indCol = grep("SLL", names(combined_ONT_WGBS))

# Save table
write.csv(combined_ONT_WGBS, "/data/SBCS-EizaguirreLab/Charley/GigaDB_Upload_Jan25/04_Methylome/WGBS_Comparison/ONT_vs_WGBS_Meth_per_CpG.csv", row.names=F)

### Run lm per individual ###

df_lm=data.frame(WGBSind=character(),Df=numeric(),F=numeric(), p=numeric())

cols = grep("SLL", colnames(combined_ONT_WGBS), value=TRUE)

for (i in cols){
  # Create and run model
  formula <- as.formula(paste0("percMeth ~ ", i, "*featureType"))
  lm <- lm(formula, data=combined_ONT_WGBS)
  results <- as.data.frame(anova(lm))
  
  # Add results to dataframe - getting 3rd row of results for interaction
  new_row <- data.frame(
    WGBSind = i,
    Df = results$Df[4],
    F = results$`F value`[3],
    p = results$`Pr(>F)`[3]
  )
  
  # Bind the new row
  df_lm <- rbind(df_lm, new_row)
}

### Run Pearson's per individual ###
dfCorWGBSONT=data.frame(WGBSind=numeric(),
                        R2=numeric(), t=numeric(), 
                        p.value=numeric(),  df=numeric())

df_exon <- dfCorWGBSONT
df_introns <- dfCorWGBSONT
df_prom <- dfCorWGBSONT
df_int <- dfCorWGBSONT

for (i in indCol){
  # Split by feature type
  exons <- combined_ONT_WGBS[combined_ONT_WGBS$featureType == "exons",]
  introns <- combined_ONT_WGBS[combined_ONT_WGBS$featureType == "introns",]
  promoters <- combined_ONT_WGBS[combined_ONT_WGBS$featureType == "promoters",]
  intergenic <- combined_ONT_WGBS[combined_ONT_WGBS$featureType == "intergenic",]
  
  # Perform correlation per feature type
  cor_exon <- cor.test(exons$percMeth, exons[ ,i])
  cor_introns <- cor.test(introns$percMeth, introns[ ,i])
  cor_prom <- cor.test(promoters$percMeth, promoters[ ,i])
  cor_int <- cor.test(intergenic$percMeth, intergenic[ ,i])
  
  # Bind into df per feature type
  df_exon = rbind(df_exon, 
                       data.frame(WGBSind=names(exons)[i],
                                  R2=round(cor_exon$estimate, 3), t=cor_exon$statistic, 
                                  p.value=cor_exon$p.value,  df=cor_exon$parameter))
  
  df_introns = rbind(df_introns, 
                  data.frame(WGBSind=names(introns)[i],
                             R2=round(cor_introns$estimate, 3), t=cor_introns$statistic, 
                             p.value=cor_introns$p.value,  df=cor_introns$parameter))
  
  df_prom = rbind(df_prom, 
                  data.frame(WGBSind=names(promoters)[i],
                             R2=round(cor_prom$estimate, 3), t=cor_prom$statistic, 
                             p.value=cor_prom$p.value,  df=cor_prom$parameter))
  
  df_int = rbind(df_int, 
                  data.frame(WGBSind=names(intergenic)[i],
                             R2=round(cor_int$estimate, 3), t=cor_int$statistic, 
                             p.value=cor_int$p.value,  df=cor_int$parameter))
}

########################################################################################################################################

#######################################
###### Per gene by feature type #######
#######################################

######## ONT vs population average WGBS methylome ########

### Prep data ###
library(dplyr)
combined_ONT_WGBS10ave_Bygene.ft <- combined_ONT_WGBS10ave %>% 
  group_by(geneInfo, featureType) %>%
  summarise(
    avg_percMeth = mean(percMeth, na.rm = TRUE),
    avg_percMeth10ave = mean(percMeth10ave, na.rm = TRUE)
  )

### Linear model ###
lm <- lm(avg_percMeth ~ avg_percMeth10ave*featureType, data=combined_ONT_WGBS10ave_Bygene.ft)
anova(lm)

### Pearson's test ###

# -> Interaction by feature type -> Pearson's correlation separately per feature type as a post-hoc test
# Subset
exons <- combined_ONT_WGBS10ave_Bygene.ft[combined_ONT_WGBS10ave_Bygene.ft$featureType == "exons",]
introns <- combined_ONT_WGBS10ave_Bygene.ft[combined_ONT_WGBS10ave_Bygene.ft$featureType == "introns",]
promoters <- combined_ONT_WGBS10ave_Bygene.ft[combined_ONT_WGBS10ave_Bygene.ft$featureType == "promoters",]
intergenic <- combined_ONT_WGBS10ave_Bygene.ft[combined_ONT_WGBS10ave_Bygene.ft$featureType == "intergenic",]

# exons
cor_exon <- cor.test(exons$avg_percMeth, exons$avg_percMeth10ave)
data.frame(t=cor_exon$statistic, df=cor_exon$parameter, 
           p.value=cor_exon$p.value, R2=round(cor_exon$estimate, 3))

# introns
cor_intron <- cor.test(introns$avg_percMeth, introns$avg_percMeth10ave)
data.frame(t=cor_intron$statistic, df=cor_intron$parameter, 
           p.value=cor_intron$p.value, R2=round(cor_intron$estimate, 3))

# promoters
cor_prom <- cor.test(promoters$avg_percMeth, promoters$avg_percMeth10ave)
data.frame(t=cor_prom$statistic, df=cor_prom$parameter, 
           p.value=cor_prom$p.value, R2=round(cor_prom$estimate, 3))

# intergenic
cor_int <- cor.test(intergenic$avg_percMeth, intergenic$avg_percMeth10ave)
data.frame(t=cor_int$statistic, df=cor_int$parameter, 
           p.value=cor_int$p.value, R2=round(cor_int$estimate, 3))

### Plot ###

scatterPlot_bygene <- ggplot(combined_ONT_WGBS10ave_Bygene.ft, 
                             aes(x=avg_percMeth, y=avg_percMeth10ave, color=featureType)) +
  geom_point(alpha=0.5, size=0.1) +
  stat_smooth(method = "lm", color="black", alpha=0.8, linewidth = 0.8) +
  labs(x="ONT: mean methylation per gene/feature type (%)", y="WGBS: mean methylation per gene/feature type (%)", color="Feature type") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  theme_classic() +
  scale_color_manual(values = c(viridis_custom_palette[5], viridis_custom_palette[1], viridis_custom_palette[2], viridis_custom_palette[3])) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16), legend.position="none") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  facet_wrap(~featureType, 
             labeller = labeller(featureType = c("exons" = "Exon", "intergenic" = "Intergenic (<10kb from TSS)", "introns" = "Intron", "promoters" = "Promoter"))) +
  theme(strip.text.x = element_text(size = 14))

ggsave(file=file.path(OUTDIR_Figures ,"All_Genes/ScatterPlot_ONT_vs_WGBS_MeanMeth_AnnotCpGs_ByGene.pdf"), 
       plot=scatterPlot_bygene , 
       width=7, height=7, units="in")

######## Individually per WGBS methylome ########

### Calculate per individual separately ###
sll_means <- combined_ONT_WGBS %>%
  group_by(geneInfo, featureType) %>%
  summarise(
    avg_percMeth = mean(percMeth, na.rm = TRUE),
    across(starts_with("SLL"), mean, na.rm = TRUE)
  )

sll_means
# colnames(sll_means)[colnames(sll_means) == 'avg_percMeth'] <- 'ONT'
write.csv(sll_means, "/data/SBCS-EizaguirreLab/Charley/GigaDB_Upload_Jan25/04_Methylome/WGBS_Comparison/ONT_vs_WGBS_Meth_per_Gene.csv", row.names=F)
saveRDS(sll_means, "/data/SBCS-EizaguirreLab/Charley/GigaDB_Upload_Jan25/04_Methylome/WGBS_Comparison/ONT_vs_WGBS_Meth_per_Gene.RDS")


### Run lm per individual ###

df_lm=data.frame(WGBSind=character(),Df=numeric(),F=numeric(), p=numeric())

for (i in cols){
  # Create and run model
  formula <- as.formula(paste0("percMeth ~ ", i, "*featureType"))
  lm <- lm(formula, data=sll_means)
  results <- as.data.frame(anova(lm))
  
  # Add results to dataframe - getting 3rd row of results for interaction
  new_row <- data.frame(
    WGBSind = i,
    Df = results$Df[4],
    F = results$`F value`[3],
    p = results$`Pr(>F)`[3]
  )
  
  # Bind the new row
  df_lm <- rbind(df_lm, new_row)
}


### Run Pearson's per individual ###

indCol = grep("SLL", names(sll_means))

dfCorWGBSONT=data.frame(WGBSind=numeric(),
                        R2=numeric(), t=numeric(), 
                        p.value=numeric(),  df=numeric())

df_exon <- dfCorWGBSONT
df_introns <- dfCorWGBSONT
df_prom <- dfCorWGBSONT
df_int <- dfCorWGBSONT

# Split by feature type
exons <- as.data.frame(sll_means[sll_means$featureType == "exons",])
introns <- as.data.frame(sll_means[sll_means$featureType == "introns",])
promoters <- as.data.frame(sll_means[sll_means$featureType == "promoters",])
intergenic <- as.data.frame(sll_means[sll_means$featureType == "intergenic",])

for (i in indCol){
  # Perform correlation per feature type
  cor_exon <- cor.test(exons$avg_percMeth, exons[ ,i])
  cor_introns <- cor.test(introns$avg_percMeth, introns[ ,i])
  cor_prom <- cor.test(promoters$avg_percMeth, promoters[ ,i])
  cor_int <- cor.test(intergenic$avg_percMeth, intergenic[ ,i])
  
  # Bind into df per feature type
  df_exon = rbind(df_exon, 
                  data.frame(WGBSind=names(exons)[i],
                             R2=round(cor_exon$estimate, 3), t=cor_exon$statistic, 
                             p.value=cor_exon$p.value,  df=cor_exon$parameter))
  
  df_introns = rbind(df_introns, 
                     data.frame(WGBSind=names(introns)[i],
                                R2=round(cor_introns$estimate, 3), t=cor_introns$statistic, 
                                p.value=cor_introns$p.value,  df=cor_introns$parameter))
  
  df_prom = rbind(df_prom, 
                  data.frame(WGBSind=names(promoters)[i],
                             R2=round(cor_prom$estimate, 3), t=cor_prom$statistic, 
                             p.value=cor_prom$p.value,  df=cor_prom$parameter))
  
  df_int = rbind(df_int, 
                 data.frame(WGBSind=names(intergenic)[i],
                            R2=round(cor_int$estimate, 3), t=cor_int$statistic, 
                            p.value=cor_int$p.value,  df=cor_int$parameter))
}
