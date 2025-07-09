#!/bin/bash -l
#SBATCH --job-name=blobtools_array     # Job name
#SBATCH --account=microalgae           # Your account
#SBATCH --partition=tier3              # Your partition
#SBATCH --time=0-01:00:00              # Time limit (D-HH:MM:SS). Adjust as needed. Blobtools create can be memory/CPU intensive.
#SBATCH --output=%x_%A_%a.log         # Standard output log (JobName_JobID_ArrayTaskID.log)
#SBATCH --error=%x_%A_%a.err          # Standard error log (JobName_JobID_ArrayTaskID.err)
#SBATCH --nodes=1                      # Run on a single node per task
#SBATCH --ntasks=1                     # One task per node
#SBATCH --cpus-per-task=8              # Number of CPU cores per task for blobtools create (adjust based on resources and file size)
#SBATCH --mem=32GB                     # Total memory for the task (adjust based on BAM file size)

## JOB ARRAY - IMPORTANT!
# Adjust this range (1-3) to match the numbers of your AKH files (e.g., 1-3 for AKH1 to AKH3)
#SBATCH --array=1-3

echo "All jobs in this array have:"
echo "- SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"
echo "- SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}"
echo "- SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN}"
echo "- SLURM_ARRAY_TASK_MAX=${SLURM_ARRAY_TASK_MAX}"
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

echo "Starting Blobtools 'create' Task: AKH${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"

# --- Configuration ---
# IMPORTANT: Adjust these paths and names to your specific setup!

# Conda environment name where blobtools is installed
CONDA_ENV_NAME="blobtools" # Make sure this matches your environment name

# Directory containing your sorted BAM files (e.g., from the previous samtools script)
BAM_DIR="/home/frvsbi/metagenomes/MH_AKH/bwa-array" # Adjust this path

# Path to your assembly FASTA file (the same one used for BWA alignment)
# This is crucial for blobtools to map reads back to contigs
ASSEMBLY_FASTA="/shared/rc/microalgae/metagenomes/assemblies/MH_AKH/MH_AKH-final.contigs.fa" # Adjust this path

# Optional: Path to a NCBI BLASTn output (or similar taxonomic classification)
# This is usually a .out or .tsv file from running blastn against a known database
# If you don't have this, you can omit the --hits argument, but results will be less informative.
BLASTN_HITS_FILE="/home/frvsbi/metagenomes/MH_AKH/MH_AKH.contigs-blastn.out" # Adjust this path if you have it

# Optional: Path to a protein file for BUSCO/Gene calling results
# This is usually a .faa or .tsv file with gene predictions
#PROTEIN_FILE="/path/to/your/protein_predictions/AKH${SLURM_ARRAY_TASK_ID}.faa" # Adjust this path if you have it

# Output directory for blobtools results for each sample
# A subdirectory will be created within this for each sample's results
BLOBTOOLS_OUTPUT_BASE_DIR="/home/frvsbi/metagenomes/MH_AKH/blobtools" # Adjust this path

# Create the base output directory if it doesn't exist
#mkdir -p "$BLOBTOOLS_OUTPUT_BASE_DIR"

# --- Conda Environment Activation ---
# This is the standard way to activate a conda environment in a script
source "$(conda info --base)/etc/profile.d/conda.sh" # Initializes conda for the shell
conda activate "$CONDA_ENV_NAME"
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate Conda environment: $CONDA_ENV_NAME. Exiting."
    exit 1
fi
echo "Conda environment '$CONDA_ENV_NAME' activated."

# --- Determine the current sample's prefix based on SLURM_ARRAY_TASK_ID ---
SAMPLE_PREFIX="AKH${SLURM_ARRAY_TASK_ID}"

INPUT_BAM="${BAM_DIR}/${SAMPLE_PREFIX}.sorted.bam"
# Each sample will have its own output directory for blobtools create
SAMPLE_OUTPUT_DIR="${BLOBTOOLS_OUTPUT_BASE_DIR}/${SAMPLE_PREFIX}"

mkdir -p "$SAMPLE_OUTPUT_DIR" # Create output directory for the current sample

echo "Processing sample: ${SAMPLE_PREFIX}"
echo "Input BAM: ${INPUT_BAM}"
echo "Output Directory: ${SAMPLE_OUTPUT_DIR}"

# --- Check if input files exist ---
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file not found: $INPUT_BAM. Skipping this task."
    exit 0 # Exit successfully for this task if file is missing
fi

if [ ! -f "$ASSEMBLY_FASTA" ]; then
    echo "Error: Assembly FASTA not found: $ASSEMBLY_FASTA. Exiting."
    exit 1
fi

# --- Run blobtools create ---
echo "Running blobtools create for ${SAMPLE_PREFIX}..."

# Construct the blobtools create command
BLOBTOOLS_CMD="blobtools create \
    --infile ${ASSEMBLY_FASTA} \
    --bam ${INPUT_BAM} \
    --out ${SAMPLE_OUTPUT_DIR}/${SAMPLE_PREFIX}"

# Add optional arguments if files exist
if [ -f "$BLASTN_HITS_FILE" ]; then
    echo "Adding --hits: ${BLASTN_HITS_FILE}"
    BLOBTOOLS_CMD+=" --hits ${BLASTN_HITS_FILE}"
fi

if [ -f "$PROTEIN_FILE" ]; then
    echo "Adding --proteins: ${PROTEIN_FILE}"
    BLOBTOOLS_CMD+=" --proteins ${PROTEIN_FILE}"
fi

# Add any other optional arguments you need, e.g., --taxonomy for a NCBI tax dump
# BLOBTOOLS_CMD+=" --taxonomy /path/to/your/ncbi_taxdump"

# Execute the command
eval "$BLOBTOOLS_CMD" # Use eval to execute the dynamically constructed command
if [ $? -ne 0 ]; then
    echo "Error: blobtools create failed for ${SAMPLE_PREFIX}. Exiting."
    conda deactivate # Deactivate Conda environment on error
    exit 1
fi
echo "blobtools create for ${SAMPLE_PREFIX} complete."

# --- Deactivate Conda Environment ---
conda deactivate
echo "Conda environment deactivated."

echo "Blobtools Task for ${SAMPLE_PREFIX} Finished!"
echo "Date: $(date)"