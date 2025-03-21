#!/bin/bash
#$ -j Y
#$ -cwd
#$ -pe smp 36
#$ -l h_rt=240:00:00
#$ -l h_vmem=10G
#$ -l highmem
#$ -m bea

#######################################################################

# Created by: Charley, 2021
# Perform de novo assembly using trimmed ONT reads with Flye 

#######################################################################

# Prep environment
module load python/3.8.5

FLYE=/data/home/btx902/privatemodules/Flye/bin/flye
DATA_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/02_Quality_Control/03_Trimmed_Reads/02_NanoFilt
OUT_DIR=/data/SBCS-EizaguirreLab/Turtle_Genome/03_De_Novo_Assembly/Flye_Assembly/asm_cov_40

# Run Flye
echo "Starting assembly: "`date`

# --asm-coverage 40 recommended for initial disjointig assembly with lower coverage
# to reduce memory consumption for large genome assemblies
python $FLYE \
  -t ${NSLOTS} \
  --nano-raw $DATA_DIR/turtle_porechop_q8minlen500.fastq.gz \
  --out-dir $OUT_DIR \
  -g 3g \
  --asm-coverage 40 \

echo "Finished assembly: "`date`
