 #!/bin/bash -l
#SBATCH --job-name=biobakery_test   # Job name for BWA alignment and Samtools processing
#SBATCH --account=microalgae             # Your account
#SBATCH --partition=tier3                # Your partition
#SBATCH --time=0-12:00:00                # Time limit (D-HH:MM:SS). Adjust as needed, BWA can be lengthy.
#SBATCH --output=%x_%A_%a.log           # Standard output log (JobName_JobID_ArrayTaskID.log)
#SBATCH --error=%x_%A_%a.err            # Standard error log (JobName_JobID_ArrayTaskID.err)
#SBATCH --ntasks=1                       # One task per node
#SBATCH --cpus-per-task=6                 # Number of CPU cores per task (BWA and Samtools can use many threads)
#SBATCH --mem=60G                # Memory per CPU. Total memory = cpus-per-task * mem-per-cpu (e.g., 16 * 8GB = 128GB)
#SBATCH --mail-user=slack:@rlk5096       # Replace with your email
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-2

# Load conda env
source /home/rlk5096/miniconda3/etc/profile.d/conda.sh
conda activate biobakery_env

# Load the sample name from file
SAMPLE_LIST="/shared/rc/biocrusts/sampleFiles/sampleTOL2_nolane.txt"
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ${SAMPLE_LIST})

# Set input paths based on sample name
FASTQ="/shared/rc/microalgae/biobakery_real/kneaddata/${SAMPLE}/${SAMPLE}.fastq.gz"
TAXON="/shared/rc/microalgae/biobakery_real/20_metaphlan/${SAMPLE}/${SAMPLE}_taxonomic_profile.tsv"
OUT_DIR="/shared/rc/microalgae/biobakery_real/humann_out/${SAMPLE}"

# Make output dir
mkdir -p "$OUT_DIR"

humann \
    --input "$FASTQ" \
    --output "$OUT_DIR" \
    --o-log /shared/rc/microalgae/biobakery_real/20_metaphlan/site2_S1.log \
    --threads 6 \
    --taxonomic-profile "$TAXON" \
    --bypass-translated-search \
    --bowtie-options "--bt2_ps very-sensitive-local" 
