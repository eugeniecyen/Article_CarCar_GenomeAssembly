#!/bin/bash
#$ -cwd
#$ -j y
#$ -m bea 
#$ -pe smp 1
#$ -l h_vmem=45G
#$ -l h_rt=240:0:0
#$ -t 1-2418

#######################################################################

# Created by: Charley, 2021
# Array per fastq file preparing raw ONT reads by trimming residual adapters 
# via PoreChop, then filtering by quality with NanoFilt

#######################################################################

### Use PoreChop to trim residual adapters ###

# Set up environment
module load gcc/8.2.0
module load python/3.8.5
source /data/SBCS-EizaguirreLab/Charley/environments/porechop_env/bin/activate

DATA_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/01_PromethION_Raw/01_Fastq_Pass
OUT_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/02_Quality_Control/04_Trimmed_Reads/01_PoreChop

# Specify file for data processing step
query_file=`sed ${SGE_TASK_ID}"q;d" turtle_rawreads_fastqc_files_list.txt`
echo $query_file

# Run PoreChop
porechop -t 4 --verbosity 1 \
-i $query_file \
-o $OUT_DIR/"porechop_${query_file:71}"

deactivate

### Use NanoFilt to filter reads by Phred >Q8 and min length 500 bp ###

# Set up environment
source /data/SBCS-EizaguirreLab/Charley/environments/nanopack_env/bin/activate

# Run NanoFilt
# NB. NanoFilt does not provide options for input or output files → use cat or redirect operators
gunzip -c $OUT_DIR/"porechop_${query_file:71}" | \
NanoFilt -l 500 -q 8 | \
gzip > $OUT_DIR/"porechop_pass_${query_file:71}"

# After array finished and QC, all fastq files were combined into a single file 
# 'turtle_porechop_q8minlen500.fastq.gz' using the cat command 
