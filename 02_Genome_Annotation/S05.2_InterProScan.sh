#!/bin/bash
#$ -pe smp 10
#$ -l h_vmem=5G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y
#$ -m bea

#######################################################################

# Created by: Charley, 2022
# Functionally characterise proteins with InterProScan

#######################################################################

module load java/17.0.0
module load python/3.10.7

INTERPROSCAN=/data/SBCS-EizaguirreLab/Charley/modules/local_modules/interproscan-5.60-92.0/interproscan.sh

$INTERPROSCAN \
-appl Gene3D,PANTHER,pfam,PIRSR,SFLD,SUPERFAMILY,TIGRFAM \
-dp -f TSV,GFF3 -goterms -iprlookup -t p \
-i CarCar_QM_v1_2021_12_Scaff_Annotation.aa.fa \
-cpu ${NSLOTS}

# With pathways
$INTERPROSCAN \
-appl Gene3D,PANTHER,pfam,PIRSR,SFLD,SUPERFAMILY,TIGRFAM \
-dp -f TSV,GFF3 -goterms -iprlookup -t p -pa \
-i CarCar_QM_v1_2021_12_Scaff_Annotation.aa.fa \
-cpu ${NSLOTS}
