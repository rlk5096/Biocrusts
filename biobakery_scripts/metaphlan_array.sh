#!/bin/bash -l
#SBATCH --job-name=metaphlan_array
#SBATCH --account=microalgae
#SBATCH --partition=tier3
#SBATCH --time=1-00:00:00
#SBATCH --output=metaphlan_%A_%a.out
#SBATCH --error=metaphlan_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --array=0-2
#SBATCH --mail-user=you@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

## Load conda environment
source /home/rlk5096/miniconda3/etc/profile.d/conda.sh
conda activate biobakery_env

## Path containing list of sample names (ex: Toledo11 for Toledo11.fastq.gz)
SAMPLE_LIST=/shared/rc/biocrusts/sampleFiles/sampleTOL1_nolane.txt


# Get sample ID from array index
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SAMPLE_LIST)

## Input is Kneaddata output zipped file
INPUT_DIR=/shared/rc/microalgae/biobakery_real/kneaddata/$SAMPLE
OUTPUT_DIR=/shared/rc/microalgae/biobakery_real/metaphlan/

# Construct file names
input_file="${INPUT_DIR}/${SAMPLE}.fastq.gz"

# Create output dir
mkdir -p ${OUTPUT_DIR}/${SAMPLE}

echo $input_file

metaphlan $input_file \
  --input_type fastq \
  --nproc 6 \
  --output_file $OUTPUT_DIR/$inputfile_taxonomic_profile.tsv \
  --bt2_ps very-sensitive-local \
  --unclassified_estimation \
  --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
  --stat_q 0.1 \
  --min_mapq_val 20 \
