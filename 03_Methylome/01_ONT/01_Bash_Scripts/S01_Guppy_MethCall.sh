#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8        # 8 cores (8 cores per GPU)
#$ -l h_rt=240:0:0   
#$ -l h_vmem=11G
#$ -l gpu=1         # request 1 GPU
#$ -l centos 		# Run before Apocrita Rocky upgrade in Feb 2025

#######################################################################

# Created by: Alice Balard, 2023 (https://github.com/alicebalard)
# Call 5mC and 5hmC base modifications from ONT reads for our reference genome

#######################################################################

# Prep environment
module load use.dev # load dev modules
module load nanopore_guppy/6.5.7 # test development version

# Specify both reference genome versions
#REF='/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc_Chr0.fasta'
REF='/data/SBCS-EizaguirreLab/Turtle_Genome/00_Final_Assemblies/CarCar_QM_v1.21.12_Sc.fasta'

# NB. To obtain the correct parameters for the config file:
# nanopore_guppy guppy_basecaller --print_workflows
# Available flowcell + kit combinations are:
# flowcell       kit               barcoding config_name                    model version
# FLO-PRO002 SQK-LSK109 dna_r9.4.1_450bps_hac_prom 2021-05-05_dna_r9.4.1_promethion_384_dd219f32

# Copy and unzip raw reads in scratch (NB. a couple were already unzipped)
raw_files_gz='/data/archive/archive-SBCS-EizaguirreLab/Turtle_Genome_Assembly/Raw_Sequencing_Data/01_PromethION_Raw/20210128_1832_1E_PAF34142_f7d42285/fast5/fast5_pass'
mkdir -p /data/scratch/btx915/raw_files
cp $raw_files_gz/* /data/scratch/btx915/raw_files/.
gunzip /data/scratch/btx915/raw_files/*.gz

input_path='/data/scratch/btx915/raw_files'
mkdir -p /data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/Guppy_outdir_2025
save_path='/data/SBCS-EizaguirreLab/Turtle_Genome/13_ONT_Methylome_Alice/output/Guppy_outdir_2025'

# Run Guppy base caller
# In config argument, modbases_5hmc_5mc_cg specifies that we call both 5hmc and 5mc in the cg context
nanopore_guppy guppy_basecaller -i $input_path -s $save_path \
-c "dna_r9.4.1_450bps_modbases_5hmc_5mc_cg_hac_prom.cfg" \
--align_ref $REF --bam_out --compress_fastq --device cuda:0


