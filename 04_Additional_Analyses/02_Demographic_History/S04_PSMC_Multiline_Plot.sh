#!/bin/bash
#$ -pe smp 4
#$ -l h_vmem=1G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y

#######################################################################

# Created by: Charley, Jul 2023
# Produce multi-line PSMC plot and output text file

#######################################################################

module load anaconda3
conda activate psmc
module load gnuplot

OUT_DIR=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/03_PSMC/02_PSMC_Output

cd $OUT_DIR

### Plot ###

# Use -R option -> also outputs text file which you can plot in R

# Macrochromosomes
psmc_plot.pl -u 1.2e-08 -g 45 -Y 8 -R \
-M "SLK063,SRR15328383" \
CarCar_macrochroms_45y_lim8_PSMC_plot \
SLK063_macrochroms.psmc \
SRR15328383_macrochroms.psmc

# Microchromosomes
psmc_plot.pl -u 1.2e-08 -g 45 -Y 8 -R \
-M "SLK063,SRR15328383" \
CarCar_microchroms_45y_lim8_PSMC_plot \
SLK063_microchroms.psmc \
SRR15328383_microchroms.psmc