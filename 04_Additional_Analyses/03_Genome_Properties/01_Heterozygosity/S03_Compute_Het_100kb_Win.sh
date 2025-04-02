#!/bin/bash
#$ -pe smp 6
#$ -l h_vmem=1G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y

#######################################################################

# Created by: Charley, Jun 2023, following workflow in Bentley et al. (2023) PNAS
# Calculate genome-wide heterozygosity for reference individual (SLK063)
# Compute heterozygosity in an array by chromosome (100kb windows, 20kb steps)

#######################################################################

module load python # 3.8.5
source /data/SBCS-EizaguirreLab/Charley/environments/numpy/bin/activate

GENOMICS_GENERAL_DIR=/data/SBCS-EizaguirreLab/Charley/scripts/genomics_general
VCF_DIR=/data/scratch/btx902/Heterozygosity/Vcfs
GENO_DIR=/data/scratch/btx902/Heterozygosity/Het_Calculation/Simon_Martin_Het/Genos
HET_DIR=/data/scratch/btx902/Heterozygosity/Het_Calculation/Simon_Martin_Het/Het_Windows

### Extract chrom IDs ###

Chrom_File=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/chrom_names.txt
CHROM=$(sed -n "${SGE_TASK_ID}p" $Chrom_File)

### Create geno file ###
python $GENOMICS_GENERAL_DIR/VCF_processing/parseVCF.py \
-i $VCF_DIR/SLK063_monomorphic.filt.vcf.gz | \
gzip > $GENO_DIR/SLK063.mono.TrimAlt.filtpass.geno.gz

# Run popgenWindows.py in het analysis mode
 python $GENOMICS_GENERAL_DIR/popgenWindows.py \
--analysis indHet \
--windType coordinate -w 100000 -s 20000 \
-g $GENO_DIR/SLK063.mono.TrimAlt.filtpass.geno.gz \
-o $HET_DIR/SLK063.het.10kb.csv.gz \
-f phased \
-T ${NSLOTS}
