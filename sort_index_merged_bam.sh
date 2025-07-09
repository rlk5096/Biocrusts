#!/bin/bash
#SBATCH --job-name=sort_index_bam
#SBATCH --output=sort_index_bam_%j.out
#SBATCH --error=sort_index_bam_%j.err
#SBATCH --partition=tier3
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=slack:                  # add your email
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END 		            # Can be BEGIN, END, FAIL
#SBATCH --mail-type=FAIL					# Optional: also notify if the job fails

# Load any required modules (if needed) - I seemed to have trouble with the spack version of samtools, however you could activate spack here or call a conda env
source /home/rlk5096/miniconda3/etc/profile.d/conda.sh
conda activate clean_samtools

# Set paths
INPUT_BAM="/shared/rc/microalgae/metagenomes/assemblies/merged.bam"
OUTPUT_SORTED="/shared/rc/microalgae/metagenomes/assemblies/bwa_samtools/merged.sorted.bam"
TEMP_DIR="/shared/rc/microalgae/metagenomes/assemblies/bwa_samtools/tmp"  # make sure this has space

# Create temp dir if not exists
mkdir -p $TEMP_DIR

# Sort BAM
echo "Sorting BAM..."
samtools sort -@ 4 -m 6G -T $TEMP_DIR/tmp_prefix -o $OUTPUT_SORTED $INPUT_BAM

# Index BAM
echo "Indexing BAM..."
samtools index $OUTPUT_SORTED

echo "Done."
