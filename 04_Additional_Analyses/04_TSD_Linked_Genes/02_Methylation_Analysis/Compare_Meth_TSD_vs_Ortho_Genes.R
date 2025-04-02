### TSD_vs_Ortho_Genes_Meth.R ###

########################################################################################################################################

# Created by Charley, 2024
# Subset and compare TSD-linked and orthogroup genes in our reference ONT methylome

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
library(viridis)

####### Set directories ######

DIR_Metadata <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/01_Working/01_Metadata"
DIR_Working <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/01_Working/02_ONT"
OUTDIR_Figures <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/02_Output/Figures/ONT"
OUTDIR_Tables <- "/data/SBCS-EizaguirreLab/Turtle_Genome/14_Epi_TSD_Genes/02_Output/Tables/ONT"
AnnoPath <- "/data/SBCS-EizaguirreLab/Turtle_Genome/10_Annotation/06_Functional_Annotation"

####### Set colour palette ######

viridis_custom_palette <- c("#0d0887", "#5b01a5", "#a82296", "#da5b69", "#fa9e3b")

########################################################################################################################################

#####################################
### Subset for genes of interest  ###
#####################################

# Subset sites associated to SD genes and orthogroups only
# From list identified and provided by James Gilbert

###### Prep and subset files ######

# Load GRanges object with gene info
myGRangeOBJ <- readRDS(file.path(DIR_Working, "GRanges_myTurtleBedSum_GeneInfoPerSite.RDS"))

# Load table of TSD-linked and orthogroup genes to subset
CarcarSeqLoc <- read.csv(file.path(DIR_Metadata, "CarcarSeqLocations.txt"), sep = "\t")

# Check for duplicates
CarcarSeqLoc %>% nrow
CarcarSeqLoc$seq_ID %>% unique %>% length

# Subset GRanges object for only genes in supplied list
sub_myGRangeOBJ = myGRangeOBJ[myGRangeOBJ$geneInfo %in% CarcarSeqLoc$seq_ID,]

# Check how many covered
sub_myGRangeOBJ$geneInfo %>% unique %>% length

# Add gene name info
sub_myGRangeOBJ$gene_name = CarcarSeqLoc[match(sub_myGRangeOBJ$geneInfo, CarcarSeqLoc$seq_ID),"gene_name"]

# Add TSD info
sub_myGRangeOBJ$isTSD <- "TSD"
sub_myGRangeOBJ$isTSD[grep("^OG*", sub_myGRangeOBJ$gene_name)] <- "nonTSD"

# Check how many covered
sub_myGRangeOBJ$gene_name[sub_myGRangeOBJ$isTSD %in% "TSD"] %>% unique %>% length
sub_myGRangeOBJ$gene_name[sub_myGRangeOBJ$isTSD %in% "nonTSD"] %>% unique %>% length

rm(myTurtleBedSum)
rm(myannotBed12)
rm(myannotGff3)
rm(CarcarSeqLoc)

########################################################################################################################################

#########################################################
### Calculate methylation across gene by feature type ###
#########################################################

sub_myGRangeOBJ <- readRDS(file = file.path(DIR_Working, "GRanges_myTurtleBedSum_GeneInfoPerSite_TSD_Ortho_Subset.RDS"))
sub_myGRangeOBJ_df <- as.data.frame(sub_myGRangeOBJ) # Convert to df

########## Calculate mean methylation and methylation count ##########

meth_per_gene_byfeatureType = sub_myGRangeOBJ_df %>% 
  group_by(gene_name, featureType, isTSD) %>%
  summarise(meanMeth = mean(fracMeth),
            nbrCpG = n(),
            nbr99CpG = sum(fracMeth>0.99),
            nbr70CpG = sum(fracMeth>0.7)) %>%
  data.frame()

########################################################################################################################################

###################################################################
### Stats: compare mean methylation between TSD genes vs orthos ###
###################################################################

# Statistically test whether SD vs non-SD genes by feature type differ in mean methylation, with a random subset of non-SD, orthogroup 
# genes matching no. of SD genes (1000 times). This must be done to not violate statistical test assumptions, since there are so many
# more orthogroups than SD genes.

########## Create function ##########

# Function that (1) formats input object for analysis, (2) runs loop to randomly subset orthogroups n times (same no. as TSD genes) and create
# statistical models from them, and (3) save p-value outputs of each iteration to a csv
# NB. Keeping ONT and WGBS functions separate, as I think it's simpler to interpret than many nested if else loops, since multiple column differences between them. 
# Options:
# myObj = df with meth data for all TSD and ortho genes to assess
# methTest = type of methylation analysis to perform 
# "meanMeth" -> perform lm stats test on mean methylation column
# "methCount70" -> perform quasipoisson GLM on meth count (CpG with >70% meth, offset all CpG)
# "methCount99" -> perform quasipoisson GLM on meth count (CpG with >99% meth, offset all CpG)
# outputFileName = name of output csv to save, as a string
# n = number of random subsamples
myMethStatsONTFun <- function(myObj, methTest , outputFileName , n){
  ###### 1. Prep dfs ######
  # Change from character to factor
  myObj$featureType <- as.factor(myObj$featureType)
  myObj$isTSD <- as.factor(myObj$isTSD)
  
  # Set up empty dataframe to bind stats outputs of each model to
  tmp_df <- as.data.frame(matrix(ncol = 4 , nrow=0))
  colnames(tmp_df) <- c("Iteration", "p_isTSD", "p_featureType", "p_interaction")
  
  ###### 2. Loop stats to run on n random subsamples of orthos ######
  for (i in 1:n){
    ### 2a. Subset random orthos ###
    # Create vector of all ortho gene names to randomly sample from
    orthos <- myObj[myObj$isTSD == "nonTSD",]
    orthos_vector <- as.vector(unique(orthos$gene_name)) 
    
    # Randomly sample of size equal to no. of TSD genes
    tsd_genes <- myObj[myObj$isTSD == "TSD",]
    num_genes <- length(unique(tsd_genes$gene_name))
    test_orthos <- sample(orthos_vector, size = num_genes) 
    length(test_orthos)
    
    # Subset df for only these orthos, and add back SD genes
    test_orthos_df <- myObj[myObj$gene_name %in% test_orthos, ] # Subset genes that match random test ortho list
    tsd_df <- myObj[myObj$isTSD == "TSD",] # Subset TSD genes only
    test_df <- rbind(tsd_df, test_orthos_df)
    rm(tsd_df)
    rm(test_orthos_df)
    
    ### 2b. Perform stats ###
    if (methTest == "meanMeth"){
      print(paste0("Running linear model stats test on mean methylation data, i=", i))
      # Set up ANOVA of lm model
      mod1 = lm(meanMeth ~ isTSD * featureType  , test_df)
      
      # Access p-values to save
      results <- as.data.frame(anova(mod1))
      p_isTSD <- results[1,5]
      p_featureType <- results[2,5]
      p_int <- results[3,5]
      
      # Convert to data frame to append to main df
      p_df <- data.frame(i, p_isTSD, p_featureType, p_int)
      
      ### (3) Bind to empty data frame
      tmp_df <- rbind(tmp_df, p_df)
      
    } else if (methTest == "methCount99"){
      print(paste0("Running QuasiPoisson GLM stats test on highly methylated site count data (high meth: >99%), i=", i))
      # Set up QuasiPoisson GLM
      poisson.model <- glm(nbr99CpG ~ isTSD * featureType + offset(log(nbrCpG)), 
                           data = test_df, family = quasipoisson(link = "log"))
      
      # Run ANOVA
      results <- as.data.frame(anova(poisson.model, test="Chisq"))
      
      # Access just p-values to save
      p_isTSD <- results[2,5]
      p_featureType <- results[3,5]
      p_int <- results[4,5]
      
      # Convert to data frame to append to main df
      p_df <- data.frame(i, p_isTSD, p_featureType, p_int)
      
      ### (3) Bind to empty data frame ###
      tmp_df <- rbind(tmp_df, p_df)
      
    } else if (methTest == "methCount70"){
      print(paste0("Running QuasiPoisson GLM stats test on highly methylated site count data (high meth: >70%), i=", i))
      # Set up QuasiPoisson GLM
      poisson.model <- glm(nbr70CpG ~ isTSD * featureType + offset(log(nbrCpG)), 
                           data = test_df, family = quasipoisson(link = "log"))
      
      # Run ANOVA
      results <- as.data.frame(anova(poisson.model, test="Chisq"))
      
      # Access just p-values to save
      p_isTSD <- results[2,5]
      p_featureType <- results[3,5]
      p_int <- results[4,5]
      
      # Convert to data frame to append to main df
      p_df <- data.frame(i, p_isTSD, p_featureType, p_int)
      
      ### (3) Bind to empty data frame ###
      tmp_df <- rbind(tmp_df, p_df)
    }
  }
  
  print(paste0("Finished ", n, " iterations. Head of output file:"))
  print(head(tmp_df))
  
  ####### 3. Save results #######
  write.csv(tmp_df, file.path(OUTDIR_Tables, outputFileName), row.names=FALSE)
  print(paste0("Saved file ", outputFileName, " to ", OUTDIR_Tables))
}


########## Run ##########

myMethStatsONTFun(meth_per_gene_byfeatureType, methTest = "meanMeth" , outputFileName = "ONT_MeanMeth_ANOVA_1000i.csv" , n=1000) # Mean meth
myMethStatsONTFun(meth_per_gene_byfeatureType, methTest = "methCount70" , outputFileName = "ONT_MethCount70_ANOVA_1000i.csv" , n=1000) # Meth count (>70%)
# myMethStatsONTFun(meth_per_gene_byfeatureType, methTest = "methCount99" , outputFileName = "ONT_MethCount99_ANOVA_1000i.csv" , n=1000) # Meth count (>99%) # NB. Not using this comparison

########## Plot distribution of results ##########

### Mean methylation ###

# Load file 
inputFile <- "ONT_MeanMeth_ANOVA_1000i.csv"
stats_df <- read.csv(file.path(OUTDIR_Tables, inputFile))

# Count how many iterations passed p=0.05
nrow(stats_df[stats_df$p_isTSD < 0.05,]) 
nrow(stats_df[stats_df$p_featureType < 0.05,])
nrow(stats_df[stats_df$p_int < 0.05,])

# Convert dataframe to long (stacked) format
stats_df_long <- reshape2::melt(stats_df, id='i',
                                variable.name = "type", 
                                value.name = "p_value")

# Designate line for p=0.05
Vertical_Line <- geom_vline(xintercept=c(0.05), linetype="dotdash", colour="black")

# Histogram plot
histPlot <- ggplot(stats_df_long, aes(p_value, color = type, fill = type)) +
  geom_histogram(bins=1000) +
  labs(x="p-value", y="Frequency") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_color_manual(values = c(viridis_custom_palette[1], viridis_custom_palette[2], viridis_custom_palette[5])) +
  scale_fill_manual(values = c(viridis_custom_palette[1], viridis_custom_palette[2], viridis_custom_palette[5])) +
  Vertical_Line +
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) + 
  facet_wrap(~type, nrow=3, scales = "free_y", 
             labeller = labeller(type = c("p_isTSD" = "Gene type (TSD or not)", "p_featureType" = "Feature type", "p_int" = "Interaction") )) +
  theme(strip.text.x = element_text(size = 12)) +
  theme(legend.position="none")

histPlot

ggsave(file=file.path(OUTDIR_Figures ,"Stats/ONT_Histogram_MeanMeth_Pvalues_1000i_Rev2.pdf"), plot=histPlot , 
       width=6.5, height=7, units="in")

### Methylation count ###

inputFile <- "ONT_MethCount70_ANOVA_1000i.csv"
stats_df <- read.csv(file.path(OUTDIR_Tables, inputFile))

# Count how many iterations passed p=0.05
nrow(stats_df[stats_df$p_isTSD < 0.05,]) 
nrow(stats_df[stats_df$p_featureType < 0.05,])
nrow(stats_df[stats_df$p_int < 0.05,])

# Convert dataframe to long (stacked) format

stats_df_long <- reshape2::melt(stats_df, id='i',
                                variable.name = "type", 
                                value.name = "p_value")

# Designate line for p=0.05
Vertical_Line <- geom_vline(xintercept=c(0.05), linetype="dotdash", colour="black")

# Histogram plot
histPlot <- ggplot(stats_df_long, aes(p_value, color = type, fill = type)) +
  geom_histogram(bins=1000) +
  labs(x="p-value", y="Frequency") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_color_manual(values = c(viridis_custom_palette[1], viridis_custom_palette[2], viridis_custom_palette[5])) +
  scale_fill_manual(values = c(viridis_custom_palette[1], viridis_custom_palette[2], viridis_custom_palette[5])) +
  Vertical_Line +
  theme_classic() + 
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) + 
  facet_wrap(~type, nrow=3, scales = "free_y", 
             labeller = labeller(type = c("p_isTSD" = "Gene type (TSD or not)", "p_featureType" = "Feature type", "p_int" = "Interaction") )) +
  theme(strip.text.x = element_text(size = 12)) +
  theme(legend.position="none")

histPlot

outputFile <- "Stats/ONT_Histogram_MethCount70_Pvalues_1000i_Rev2.pdf"

ggsave(file=file.path(OUTDIR_Figures , outputFile), plot=histPlot , 
       width=6.5, height=7, units="in")


########################################################################################################################################

######################
### Visualise data ###
######################

plot_df <- meth_per_gene_byfeatureType
plot_df$meanMeth <- plot_df$meanMeth*100 # Change to % so it matches WGBS

### Violin plot of mean methylation ###
violinPlot <- ggplot(plot_df, aes(y=meanMeth, fill=isTSD, x=featureType)) +
  geom_violin(alpha=0.8) +
  labs(x="Feature type", y="Mean methylation (%)", fill="Gene type") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +
  scale_x_discrete(labels=c('Exon', 'Intergenic', 'Intron', "Promoter")) +
  scale_fill_manual(labels = c("non-TSD", "TSD"), values = c(viridis_custom_palette[1], viridis_custom_palette[5])) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

violinPlot

ggsave(file=file.path(OUTDIR_Figures , "By_Gene_AllOrthos/ONT_MeanMeth_ViolinPlot.pdf"), plot=violinPlot, width=8, height=8, units="in")

##### Violin plot of highly methylated sites proportion #####

violinPlot <- ggplot(plot_df, aes(y=nbr70CpG/nbrCpG, fill=isTSD, x=featureType)) +
  geom_violin(alpha=0.8) +
  labs(x="Feature type", y="Ratio of highly methylated (>70%) sites", fill="Gene type") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_x_discrete(labels=c('Exon', 'Intergenic', 'Intron', "Promoter")) +
  scale_fill_manual(labels = c("non-TSD", "TSD"), values = c(viridis_custom_palette[1], viridis_custom_palette[5])) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

violinPlot

ggsave(file=file.path(OUTDIR_Figures , "By_Gene_AllOrthos/ONT_MethCount70_ViolinPlot.pdf"), plot=violinPlot, width=8, height=8, units="in")

### Density plot for just 1 random subset of orthos ###

# Density plot affected by total numbers, so looks different between orthogroups and TSD genes
# Plot for just a random subset of orthos
meth_per_gene_byfeatureType$featureType <- as.factor(meth_per_gene_byfeatureType$featureType)
meth_per_gene_byfeatureType$isTSD <- as.factor(meth_per_gene_byfeatureType$isTSD)

# Create vector of all ortho gene names to randomly sample from
orthos <- meth_per_gene_byfeatureType[meth_per_gene_byfeatureType$isTSD == "nonTSD",]
orthos_vector <- as.vector(unique(orthos$gene_name))
# Randomly sample of size equal to no. of TSD genes
tsd_genes <- meth_per_gene_byfeatureType[meth_per_gene_byfeatureType$isTSD == "TSD",]
num_genes <- length(unique(tsd_genes$gene_name))
test_orthos <- sample(orthos_vector, size = num_genes)
length(test_orthos)
# Subset df for only these orthos, and add back SD genes
test_orthos_df <- meth_per_gene_byfeatureType[meth_per_gene_byfeatureType$gene_name %in% test_orthos, ] # Subset genes that match random test ortho list
tsd_df <- meth_per_gene_byfeatureType[meth_per_gene_byfeatureType$isTSD == "TSD",] # Subset TSD genes only
test_df <- rbind(tsd_df, test_orthos_df)
rm(tsd_df)
rm(test_orthos_df)

# Density plot: mean methylation
densityPlot <- ggplot(test_df, aes(meanMeth, fill = isTSD)) +
  geom_density(alpha=0.7, color = "black") +
  labs(x="Mean methylation (%)", y="Density", fill="Gene type") +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_fill_manual(labels = c("non-TSD", "TSD"), values = c(viridis_custom_palette[5], viridis_custom_palette[1])) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  facet_wrap(~featureType, scales="free_y",
             labeller = labeller(featureType = c("exons" = "Exon", "intergenic" = "Intergenic", "introns" = "Intron", "promoters" = "Promoter"))) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

densityPlot

ggsave(file=file.path(OUTDIR_Figures , "By_Gene_AllOrthos/ONT_MeanMeth_DensityPlot_ByFeatureType_OrthoSubset.pdf"), 
       plot=densityPlot, width=8, height=7, units="in")

# Density plot: methylation count
densityPlot <- ggplot(test_df, aes(nbr70CpG/nbrCpG, fill = isTSD)) +
  geom_density(alpha=0.7, color = "black") +
  labs(x="Ratio of highly methylated (>70%) sites", y="Density", fill="Gene type") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(labels = c("non-TSD", "TSD"), values = c(viridis_custom_palette[5], viridis_custom_palette[1])) +
  theme_classic() +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14), 
        legend.text = element_text(size=14), legend.title = element_text(size=16)) + 
  facet_wrap(~featureType, scales="free_y",
             labeller = labeller(featureType = c("exons" = "Exon", "intergenic" = "Intergenic", "introns" = "Intron", "promoters" = "Promoter"))) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

densityPlot

ggsave(file=file.path(OUTDIR_Figures , "By_Gene_AllOrthos/ONT_MethCount70_DensityPlot_ByFeatureType_OrthoSubset.pdf"), 
       plot=densityPlot, width=8, height=7, units="in")

