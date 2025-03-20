#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -m bea  
#$ -l h_rt=240:0:0
#$ -l h_vmem=45G

#######################################################################

# Created by: Charley, 2021
# Prepare raw ONT reads by trimming residual adapters using PoreChop,
# then filter by quality using NanoFilt

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

DATA_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/02_Quality_Control/04_Trimmed_Reads/01_PoreChop
OUT_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/02_Quality_Control/04_Trimmed_Reads/02_NanoFilt

# Run NanoFilt
# NB. NanoFilt does not provide options for input or output files → use cat or redirect operators
gunzip -c $DATA_DIR/combined_porechop_PAF34142_pass_a5072ac4.fastq.gz | \
NanoFilt -l 500 -q 8 | \
gzip > $OUT_DIR/nanofilt_q8minlen500_porechop.fastq.gz
