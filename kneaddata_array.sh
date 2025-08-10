#!/bin/bash -l
#SBATCH --job-name=kneaddata_array
#SBATCH --account=microalgae
#SBATCH --partition=tier3
#SBATCH --time=0-04:00:00
#SBATCH --output=kneaddata_%A_%a.out
#SBATCH --error=kneaddata_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --array=0-2
#SBATCH --mail-user=you@example.com
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda env
source /home/rlk5096/miniconda3/etc/profile.d/conda.sh
conda activate biobakery_env

# Set paths
DATABASE=/shared/rc/biocrusts/biobakery_workflows_databases/kneaddata/genome_kneaddata
INPUT_DIR=/shared/rc/biocrusts/Ohio_Samples/raw_reads
OUTPUT_DIR=/shared/rc/microalgae/biobakery_real/kneaddata
SAMPLE_LIST=/shared/rc/biocrusts/sampleFiles/sampleTOL1_nolane.txt

# Sample ID from array index
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SAMPLE_LIST)

# Construct file names
R1="${INPUT_DIR}/${SAMPLE}_R1.fastq"
R2="${INPUT_DIR}/${SAMPLE}_R2.fastq"

# Make output dir
mkdir -p ${OUTPUT_DIR}/${SAMPLE}

echo $R1
echo $R2

kneaddata \
  --input1 $INPUT_DIR/$R1 \
  --input2 $INPUT_DIR/$R2 \
  --output $OUTPUT_DIR/$SAMPLE \
  --reference-db $DATABASE \
  --threads 8 \
  --remove-intermediate-output \
  --cat-final-output

#renaming and zipping combined final output
cp ${SAMPLE}_R1_kneaddata.fastq ${SAMPLE}.fastq
gzip ${SAMPLE}.fastq
