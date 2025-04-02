#!/bin/bash
#$ -pe smp 6
#$ -l h_vmem=1G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, Jun 2023, following workflow in Bentley et al. (2023) PNAS
# Calculate genome-wide heterozygosity for reference individual (SLK063)
# Combine and filter VCFs

#######################################################################

module load java/11.0.2
module load htslib

VCF_DIR=/data/scratch/btx902/Heterozygosity/Vcfs
REPORTS_DIR=/data/scratch/btx902/Heterozygosity/Vcfs/All_Active_Sites/Reports
ASSEMBLY=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc_Chr0.fasta

cd $VCF_DIR

### Merge VCFs ###

PICARD=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/picard/build/libs/picard.jar

# Merge VCFs
java -jar $PICARD GatherVcfs \
    I=SLK063_ragtag_chr1.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr2.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr3.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr4.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr5.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr6.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr7.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr8.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr9.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr10.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr11.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr12.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr13.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr14.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr15.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr16.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr17.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr18.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr19.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr20.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr21.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr22.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr23.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr24.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr25.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr26.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr27.mono.TrimAlt.vcf.gz \
    I=SLK063_ragtag_chr28.mono.TrimAlt.vcf.gz \    
O=SLK063.mono.TrimAlt.vcf

# bgzip a vcf and create tabix index
bgzip -@ ${NSLOTS} SLK063.mono.TrimAlt.vcf && \
tabix -p vcf SLK063.mono.TrimAlt.vcf.gz

### Hard filter the raw merged VCF ###

# NB. Mean depth of sample was calculated to determine depth filters using: 
# vcftools --gzvcf SLK063.mono.TrimAlt.vcf.gz \
# --depth --out MeanDepth_mono_TrimAlt

module load gatk # v.4.2.6.1
module load vcftools

# Generate VCF with filter flags added
# Cutoff values in JEXL filter expressions must be given as doubles (0.0 not 0)
gatk VariantFiltration \
-R $ASSEMBLY \
-V SLK063.mono.TrimAlt.vcf.gz \
-O SLK063.mono.TrimAlt.flag.vcf.gz \
  --filter-expression "MQ < 50.0" --filter-name "badMQ" \
  --filter-expression "QD < 2.0" --filter-name "badQD" \
  --filter-expression "FS > 60.0" --filter-name "badFS" \
  --filter-expression "SOR > 3.0" --filter-name "badSOR" \
  --filter-expression "MQRankSum < -12.5" --filter-name "badMQRankSum" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "badReadPosRankSum" \
  --filter-expression "DP < 18.7" --filter-name "lowCov" \
  --filter-expression "DP > 112.2" --filter-name "highCov"

# Subset sites that pass all filters

gatk SelectVariants \
-R $ASSEMBLY \
-V SLK063.mono.TrimAlt.flag.vcf.gz \
--exclude-filtered true \
-O SLK063_monomorphic.filt.vcf.gz

### Generate filtering stats ###

echo -e "\n### Generate vcftools filtering summary ###\n"

vcftools --gzvcf SLK063_monomorphic.filt.vcf.gz --FILTER-summary \
--out $REPORTS_DIR/GATKFiltSummary.mono.TrimAlt

echo -e "\n### Generate no. of snps file ###\n"

zgrep -v "#" SLK063_monomorphic.filt.vcf.gz \
| wc -l > $REPORTS_DIR/NoOfSnps.mono.filtpass.txt


echo -e "\n### Check proportion of missing data ###\n"

vcftools --gzvcf SLK063_monomorphic.filt.vcf.gz \
--missing-indv \
--out $REPORTS_DIR/Missingness.mono.filtpass









