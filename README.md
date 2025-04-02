## Summary

This repository contains scripts associated with the research article: 

**Chromosome-level genome assembly and methylome profile yield insights for the conservation of endangered loggerhead sea turtles**

Authors: Eugenie C. Yen, James D. Gilbert, Alice Balard, Albert Taxonera, Kirsten Fairweather, Heather L. Ford, Doko-Miles J. Thorburn, Stephen J. Rossiter, José M. Martín-Durán, Christophe Eizaguirre

Raw sequencing reads and genome assemblies are available under study accession PRJEB79015
<br/><br/>
## Contents

### File summary

Scripts are split into 4 main directories: 

#### 01_Genome_Assembly
Scripts used for genome assembly
* `S01_ONT_Read_Trimming.sh`: run array per FASTQ file trimming raw ONT reads for residual adapters and quality
* `S02_Flye_De_Novo_Assembly.sh`: perform de novo assembly of trimmed ONT reads with Flye
* `S03_Polish_Medaka.sh`: create polished consensus sequence with Medaka
* `S04.1_Trim_Illumina_Reads.sh`: trim raw Illumina reads for residual adapters and quality
* `S04.2_Prep_Bam_Alignments_Round1.sh`: prepare analysis-ready BAM files from trimmed Illumina reads for polishing with Pilon (round 1)
* `S04.3_Polish_Pilon_Round1.sh`: polish assembly with Pilon (round 1)
* `S04.4_Prep_Bam_Alignments_Round2.sh`: prepare analysis-ready BAM files for polishing with Pilon (round 2)
* `S04.5_Polish_Pilon_Round2.sh`: polish assembly with Pilon (round 2)
* `S05_Purge_Dups.sh`: haploidise genome assembly with Purge_Dups
* `S06_Ref_Guided_Scaffolding.sh`: reference-guided scaffolding of contigs against the 'GSC_CCare_1.0' loggerhead assembly (Chang et al., 2023) with RagTag
* `S07_Assemble_Mito_Genome.sh`: extract, assemble and annotate the mitochondrial genome with MitoZ
<br/><br/>
#### 02_Genome_Annotation
Scripts used for repeat identification and genome annotation
* `S01.1_RepeatModeler.sh`: run RepeatModeler for our reference assembly
* `S01.2_Filter_Bona_Fide_Genes.sh`: remove potential bona fide genes from RepeatModeler library using the curated proteome for 'rCheMyd1.pri.v2'  (Bentley et al., 2023)
* `S01.3_Rename_TEclass.sh`: rename repeats with TEclass categories
* `S01.4_RepeatMasker.sh`: run RepeatMasker with custom repeat library to produce soft-masked genome
* `S02.1_Align_Transcriptomes_GMAP.sh`: align transcriptome hints to reference assembly with GMAP
* `S02.2a_Generate_STAR_Index.sh`: generate STAR index for reference assembly
* `S02.2b_Align_RNAseq_STAR_Array.sh`: run array aligning all RNAseq reads to genome with STAR 
* `S02.3a_Curate_Intron_Junctions_Array.sh`: run array per BAM file to curate intron junctions with Portcullis
* `S02.3b_Merge_Intron_Junctions.sh`: merge curated intron junction bed files with junctools
* `S02.4a_Mikado_Configure.sh`: prepare Mikado configuration file
* `S02.4b_Mikado_Prepare.sh`: prepare sorted, non-redundant GTF with all of the input assemblies
* `S02.4c_Generate_Homology_Info.sh`: generate homology info for Mikado predicted transcripts using the SwissProt database
* `S02.4d_Calculate_ORFs.sh`: predict open reading frames and define coding regions for each transcript with Transdecoder
* `S02.4e_Mikado_Consensus.sh`: take all evidences to create consensus set of best-fit gene models with Mikado
* `S03.1_Run_BRAKER1.sh`: perform gene prediction using RNA-Seq hints with BRAKER1
* `S04.1_PASA_Load_Mikado_Loci.sh`: load Mikado loci into a PASA database
* `S04.2_PASA_Load_Augustus_Loci.sh`: update PASA database with Augustus loci
* `S04.3_PASA_Annot_Compare_Round1.sh`: round 1 of PASA annotation comparison and updates to incorporate transcript alignnments into gene structures
* `S04.4_PASA_Annot_Compare_Round2.sh`: round 2 of PASA annotation comparison and updates to incorporate transcript alignnments into gene structures
* `S04.5_PASA_Annot_Compare_Round3.sh`: round 3 of PASA annotation comparison and updates to incorporate transcript alignnments into gene structures
* `S04.6_Clean_Annot.sh`: clean up annotation with AGAT
* `S04.7_Rename_ENSEMBL.sh`: rename genes in ENSEMBL format with prefix "CarCarScaff"
* `S05.1_SwissProt_Prot_Homology.sh`: assign functional annotation to genes via homology against the SwissProt database
* `S05.2_InterProScan.sh`: add GO terms with InterProScan
<br/><br/>
#### 03_Methylome
**ONT**

Scripts used to produce the reference ONT methylome

01_Bash_Scripts:
* `S01_Guppy_MethCall.sh`: call 5mC and 5hmC base modifications from ONT reads for the reference assembly with Guppy
* `S02.1_Merge_Guppy_Bams.sh`: merge methylation call BAM files
* `S02.2_Output_Bam_Stats.sh`: count reads and output stats for BAM file
* `S03_Modkit_BamtoBed_Destrand.sh`: convert methylation call to BED file and destrand across CpG sites

02_R_Scripts:
* `R01_methylStat_methylKit.R`: generate ONT methylation call stats and merge 5mC/5hmC sites, since WGBS cannot distinguish between them
* `R02_Annotate_Gene_Assoc_ONT.R`: annotate and extract all CpGs associated with genes
<br/><br/>

**WGBS**

Scripts used to process WGBS data for ten additional individuals and comparisons against the reference ONT methylome

01_Bash_Scripts:
* `S01.1_PrepObjects_ChromArray_75pc_Reloc3_Adults.sh`: bash submission script for R01.1_PrepObjects_ChromArray_75pc_Reloc3_Adults.R
* `S01.2_MergeChrom_75pc_Reloc3_Adults.sh`: bash submission script for R01.2_MergeChrom_75pc_Reloc3_Adults.R

02_R_Scripts: 
* `R01.1_PrepObjects_ChromArray_75pc_Reloc3_Adults.R`: run array by chromosome to create methylBase objects for CpGs covered in >75% of 10 WGBS methylomes
* `R01.2_MergeChrom_75pc_Reloc3_Adults.R`: merge methylBase files by chromosome into a whole genome file
* `R02_Annotate_Gene_Assoc_WGBS.R`: annotate and extract all CpGs associated with genes
* `R03_Compare_ONT_WGBS_Methylomes.R`: compare methyation per CpG and per gene by feature type between the ONT reference methylome 10 additional WGBS methylomes
<br/><br/>

#### 04_Additional_Analyses

**01_Whole_Genome_Alignment**

Scripts used to generate whole genome alignments for analysing genome-wide synteny between sea turtle species
* `Sea_Turtle_Genome_Alignment.sh`: produce minimap2 alignments of our reference assembly and chromosome-level assemblies for other sea turtle species for input into online DGENIES portal
<br/><br/>

**02_Demographic_History**

Scripts used to perform PSMC analysis of our reference individual SLK063 and Brazilian individual SAMN20502673 (Vilaça et al., 2021)
* `S01a_Prep_Bam_SLK063.sh`: prep analysis-ready BAM files from raw Illumina reads for reference individual SLK063
* `S01b_Prep_Bam_SAMN20502673.sh`: prep analysis-ready BAM files from raw Illumina reads for individual SAMN20502673
* `S02a_Prep_Consensus_Seq_SLK063_ArrayChrom.sh`: call consensus sequence in an array by chromosome for SLK063
* `S02b_Prep_Consensus_Seq_SAMN20502673_ArrayChrom.sh`: call consensus sequence in an array by chromosome for SAMN20502673
* `S03_Run_PSMC_SLK063.sh`: run PSMC for SLK063
* `S03_Run_PSMC_SAMN20502673.sh`: run PSMC for SAMN20502673
* `S04_PSMC_Multiline_Plot.sh`: output multi-line PSMC plot and text file
<br/><br/>

**03_Genome_Properties**

Scripts used to compute heterozygosity and % GC content for our reference assembly

01_Heterozygosity:
* `S01_Compute_Het_GATK_ArrayChrom.sh`: run an array by chromosome performing variant calling with monomorphic sites with GATK
* `S02_Compute_Het_Combine_VCFs.sh`: combine chromosome VCFs into a single gneome and filter by quality
* `S03_Compute_Het_100kb_Win.sh`: compute heterozygosity in 100kb windows

02_GC_Content: 
* `Compute_GC_Content_100kb_Win.sh`: compute GC content in 100kb windows
<br/><br/>

**04_TSD_Linked_Genes**

01_Synteny_Analysis:

Scripts used to identify TSD-linked genes in our reference assembly and the hawksbill assembly (Guo et al. 2023), based on list compiled by Bentley et al. (2023)
* `S01_Identify_TSD_Gene_Locations.sh`: identify TSD-linked genes in our reference assembly and the hawksbill assembly
* `R02.1_Check_TSD_Gene_Locations_CC.R`: check TSD-linked genes identified in our reference assembly
* `R02.2_Check_TSD_Gene_Locations_EI.R`: check TSD-linked genes identified in the hawksbill assembly
* `R03_Plot_TSD_Gene_Locations.R`: plot the chromosomal positions of TSD-linked genes between sea turtle species with circos

02_Methylation_Analysis:

Scripts used to compare methylation 
* `S01_Run_Orthofinder.sh`: identify single-copy homologues between the loggerhead, green, leatherback and hawksbill sea turtle species with OrthoFinder
* `Compare_Meth_TSD_vs_Ortho_Genes.R`: compare methylation across gene feature types between TSD-linked versus orthogroup genes in the reference ONT methylome

03_StringDB_Analysis:
* `Extract_TSD_Proteins_For_StringDB.sh`: extract TSD-linked proteins for input into online StringDB portal
<br/><br/>
### File structure overview
```
.
├── 01_Genome_Assembly
│   ├── S01_ONT_Read_Trimming.sh
│   ├── S02_Flye_De_Novo_Assembly.sh
│   ├── S03_Polish_Medaka.sh
│   ├── S04.1_Trim_Illumina_Reads.sh
│   ├── S04.2_Prep_Bam_Alignments_Round1.sh
│   ├── S04.3_Polish_Pilon_Round1.sh
│   ├── S04.4_Prep_Bam_Alignments_Round2.sh
│   ├── S04.5_Polish_Pilon_Round2.sh
│   ├── S05_Purge_Dups.sh
│   ├── S06_Ref_Guided_Scaffolding.sh
│   └── S07_Assemble_Mito_Genome.sh
├── 02_Genome_Annotation
│   ├── S01.1_RepeatModeler.sh
│   ├── S01.2_Filter_Bona_Fide_Genes.sh
│   ├── S01.3_Rename_TEclass.sh
│   ├── S01.4_RepeatMasker.sh
│   ├── S02.1_Align_Transcriptomes_GMAP.sh
│   ├── S02.2a_Generate_STAR_Index.sh
│   ├── S02.2b_Align_RNAseq_STAR_Array.sh
│   ├── S02.3a_Curate_Intron_Junctions_Array.sh
│   ├── S02.3b_Merge_Intron_Junctions.sh
│   ├── S02.4a_Mikado_Configure.sh
│   ├── S02.4b_Mikado_Prepare.sh
│   ├── S02.4c_Generate_Homology_Info.sh
│   ├── S02.4d_Calculate_ORFs.sh
│   ├── S02.4e_Mikado_Consensus.sh
│   ├── S03.1_Run_BRAKER1.sh
│   ├── S04.1_PASA_Load_Mikado_Loci.sh
│   ├── S04.2_PASA_Load_Augustus_Loci.sh
│   ├── S04.3_PASA_Annot_Compare_Round1.sh
│   ├── S04.4_PASA_Annot_Compare_Round2.sh
│   ├── S04.5_PASA_Annot_Compare_Round3.sh
│   ├── S04.6_Clean_Annot.sh
│   ├── S04.7_Rename_ENSEMBL.sh
│   ├── S05.1_SwissProt_Prot_Homology.sh
│   └── S05.2_InterProScan.sh
├── 03_Methylome
│   ├── 01_ONT
│   │   ├── 01_Bash_Scripts
│   │   │   ├── S01_Guppy_MethCall.sh
│   │   │   ├── S02.1_Merge_Guppy_Bams.sh
│   │   │   ├── S02.2_Output_Bam_Stats.sh
│   │   │   └── S03_Modkit_BamtoBed_Destrand.sh
│   │   └── 02_R_Scripts
│   │       ├── R01_methylStat_methylKit.R
│   │       └── R02_Annotate_Gene_Assoc_ONT.R
│   └── 02_WGBS
│       ├── 01_Bash_Scripts
│       │   ├── S01.1_PrepObjects_ChromArray_75pc_Reloc3_Adults.sh
│       │   └── S01.2_MergeChrom_75pc_Reloc3_Adults.sh
│       └── 02_R_Scripts
│           ├── R01.1_PrepObjects_ChromArray_75pc_Reloc3_Adults.R
│           ├── R01.2_MergeChrom_75pc_Reloc3_Adults.R
│           ├── R02_Annotate_Gene_Assoc_WGBS.R
│           └── R03_Compare_ONT_WGBS_Methylomes.R
├── 04_Additional_Analyses
│   ├── 01_Whole_Genome_Alignment
│   │   └── Sea_Turtle_Genome_Alignment.sh
│   ├── 02_Demographic_History
│   │   ├── S01a_Prep_Bam_SLK063.sh
│   │   ├── S01b_Prep_Bam_SAMN20502673.sh
│   │   ├── S02a_Prep_Consensus_Seq_SLK063_ArrayChrom.sh
│   │   ├── S02b_Prep_Consensus_Seq_SAMN20502673_ArrayChrom.sh
│   │   ├── S03_Run_PSMC_SAMN20502673.sh
│   │   ├── S03_Run_PSMC_SLK063.sh
│   │   └── S04_PSMC_Multiline_Plot.sh
│   ├── 03_Genome_Properties
│   │   ├── 01_Heterozygosity
│   │   │   ├── S01_Compute_Het_GATK_ArrayChrom.sh
│   │   │   ├── S02_Compute_Het_Combine_VCFs.sh
│   │   │   └── S03_Compute_Het_100kb_Win.sh
│   │   └── 02_GC_Content
│   │       └── Compute_GC_Content_100kb_Win.sh
│   └── 04_TSD_Linked_Genes
│       ├── 01_Synteny_Analysis
│       │   ├── R02.1_Check_TSD_Gene_Locations_CC.R
│       │   ├── R02.2_Check_TSD_Gene_Locations_EI.R
│       │   ├── R03_Plot_TSD_Gene_Locations.R
│       │   └── S01_Identify_TSD_Gene_Locations.sh
│       ├── 02_Methylation_Analysis
│       │   ├── Compare_Meth_TSD_vs_Ortho_Genes.R
│       │   └── S01_Run_Orthofinder.sh
│       └── 03_StringDB_Analysis
│           └── Extract_TSD_Proteins_For_StringDB.sh
└── README.md
```
