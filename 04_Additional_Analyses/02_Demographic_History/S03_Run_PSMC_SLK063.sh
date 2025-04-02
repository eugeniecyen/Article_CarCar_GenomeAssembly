#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=4G
#$ -l h_rt=240:00:00
#$ -cwd
#$ -j y

#######################################################################

# Created by: Charley, Jul 2023
# Run PSMC
# For reference individual (SLK063)

#######################################################################

DATA_DIR=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/03_PSMC/01_Consensus_Sequences
OUT_DIR=/data/SBCS-EizaguirreLab/Charley/01_Genome_Paper_Analyses/03_PSMC/02_PSMC_Output
ASSEMBLY=/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc_Chr0.fasta

### Merge chromosome files into a single consensus sequence ###

cd $DATA_DIR

# NB. Pre-separated chromosome files into macrochromosomes and microchromosomes

# For macrochromosomes
cat macrochroms/SLK063_ragtag_chr*.fq > macrochroms/SLK063_macrochroms_consensus.fq
gzip macrochroms/SLK063_macrochroms_consensus.fq
mv macrochroms/SLK063_macrochroms_consensus.fq.gz $OUT_DIR

# For microchromosomes
cat microchroms/SLK063_ragtag_chr*.fq > microchroms/SLK063_microchroms_consensus.fq
gzip microchroms/SLK063_microchroms_consensus.fq
mv microchroms/SLK063_microchroms_consensus.fq.gz $OUT_DIR

### Convert FASTQ file to input format for PSMC ###

module load anaconda3
conda activate psmc

# Convert FASTQ files
fq2psmcfa \
SLK063_macrochroms_consensus.fq.gz > SLK063_macrochroms_consensus.psmcfa
 
fq2psmcfa \
SLK063_microchroms_consensus.fq.gz > SLK063_microchroms_consensus.psmcfa

# Split for bootstrapping
splitfa SLK063_macrochroms_consensus.psmcfa \
> SLK063_macrochroms_consensus_split.psmcfa

splitfa SLK063_microchroms_consensus.psmcfa \
> SLK063_microchroms_consensus_split.psmcfa

### Run PSMC ###

module load gnuplot

cd $OUT_DIR

# Create function that runs PSMC for a chromosome type (macro or micro)
run_psmc_analysis() {
    local chrom_type=$1
    local output_prefix="SLK063_${chrom_type}chroms"
    local consensus_file="$DATA_DIR/SLK063_${chrom_type}chroms_consensus.psmcfa"
    local split_file="$DATA_DIR/SLK063_${chrom_type}chroms_consensus_split.psmcfa"
    
    # Run PSMC
    psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o "${output_prefix}.psmc" "$consensus_file"
    
    # Run 100 bootstraps
    seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" \
    -o round-{}.psmc "$split_file" | sh
    
    # Merge run with bootstraps
    cat "${output_prefix}.psmc" round-*.psmc > "${output_prefix}_combined_100bs.psmc"
    
    # Plot with specified generation time and mutation rate
    psmc_plot.pl -u 1.2e-08 -g 45 \
    "${output_prefix}_combined_100bs_45y_plot" \
    "${output_prefix}_combined_100bs.psmc"
}

# Run analysis for both chromosome types
for chrom_type in "macro" "micro"; do
    run_psmc_analysis "$chrom_type"
done




