#!/bin/bash -l
#SBATCH --job-name=BWA-Samtools-array    # Job name for BWA alignment and Samtools processing
#SBATCH --account=microalgae             # Your account
#SBATCH --partition=tier3                # Your partition
#SBATCH --time=2-00:00:00                # Time limit (D-HH:MM:SS). Adjust as needed, BWA can be lengthy.
#SBATCH --output=%x_%A_%a.log           # Standard output log (JobName_JobID_ArrayTaskID.log)
#SBATCH --error=%x_%A_%a.err            # Standard error log (JobName_JobID_ArrayTaskID.err)
#SBATCH --nodes=1                        # Run on a single node per task
#SBATCH --ntasks=1                       # One task per node
#SBATCH --cpus-per-task=16               # Number of CPU cores per task (BWA and Samtools can use many threads)
#SBATCH --mem-per-cpu=8g                 # Memory per CPU. Total memory = cpus-per-task * mem-per-cpu (e.g., 16 * 8GB = 128GB)
#SBATCH --mail-user=slack:@rlk5096       # Replace with your email
#SBATCH --mail-type=BEGIN,END,FAIL

## JOB ARRAY - IMPORTANT!
# Adjust this range (1-3) to match the numbers of your AKH files (e.g., 1-91 for AKH1 to AKH91)
#SBATCH --array=1-3 # Example for AKH1, AKH2, AKH3. Adjust this range to your actual sample numbers.

echo "All jobs in this array have:"
echo "- SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"
echo "- SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}"
echo "- SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN}"
echo "- SLURM_ARRAY_TASK_MAX=${SLURM_ARRAY_TASK_MAX}"
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"

# echo "Starting BWA Alignment and Samtools Processing Task: MH_TOL${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"

# --- Configuration ---
# IMPORTANT: Adjust these paths and names to your specific setup!

# Define the FULL PATH to your samtools executable
SAMTOOLS_EXE="/.autofs/tools/spack/opt/spack/linux-rhel9-skylake_avx512/gcc-12.3.1/samtools-1.19.2-op73g5maabu3h7plhb2aq74aa2pfuodx/bin/samtools"

# Path to your reference genome FASTA file
REFERENCE_GENOME="/shared/rc/microalgae/metagenomes/assemblies/MH_TOL5/final.contigs.fa"

# Directory containing your FASTQ.gz files
FASTQ_DIR="/shared/rc/biocrusts/Ohio_Samples/trimmed_reads" # Adjust this path

# Output directory for BAM files
OUTPUT_DIR="/shared/rc/microalgae/metagenomes/assemblies/MH_TOL5/bwa_samtools" # Adjust this path

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# --- Load Modules ---
# Ensure your environment is set up. This might activate BWA and other tools.
spack env activate default-genomics-x86_64-25032001

# --- Determine the current sample's prefix based on SLURM_ARRAY_TASK_ID ---
SAMPLE_PREFIX="Toledo5${SLURM_ARRAY_TASK_ID}"

R1_FILE="${FASTQ_DIR}/${SAMPLE_PREFIX}_R1.fastq.gz"
R2_FILE="${FASTQ_DIR}/${SAMPLE_PREFIX}_R2.fastq.gz"
# OUTPUT_SAM="${OUTPUT_DIR}/${SAMPLE_PREFIX}.sam" # No longer explicitly creating SAM file
OUTPUT_UNSORTED_BAM="${OUTPUT_DIR}/${SAMPLE_PREFIX}.unsorted.bam" # Intermediate unsorted BAM
SORTED_BAM="${OUTPUT_DIR}/${SAMPLE_PREFIX}.sorted.bam" # Final sorted BAM

# --- BWA Index Reference Genome (ONE-TIME SETUP) ---
# This step only needs to be run once for your reference genome.
# It's good practice to include it in the script, but if you're certain it's indexed,
# you can comment this block out. It will only run if the .bwt file is missing.
echo "Checking/Creating BWA index for reference genome..."
if [ ! -f "${REFERENCE_GENOME}.bwt" ]; then
    echo "BWA index not found. Creating index for ${REFERENCE_GENOME}..."
    bwa index "$REFERENCE_GENOME"
    if [ $? -ne 0 ]; then
        echo "Error: BWA index creation failed. Exiting."
        exit 1
    fi
    echo "BWA index created."
else
    echo "BWA index found. Skipping indexing."
fi

# --- Check if input FASTQ files exist ---
if [ ! -f "$R1_FILE" ]; then
    echo "Error: R1 file not found: $R1_FILE. Exiting for this task."
    exit 1
fi
if [ ! -f "$R2_FILE" ]; then
    echo "Error: R2 file not found: $R2_FILE. Exiting for this task."
    exit 1
fi

echo "Processing sample: ${SAMPLE_PREFIX}"
echo "R1: ${R1_FILE}"
echo "R2: ${R2_FILE}"

# --- 1. BWA Alignment and direct SAM to BAM conversion ---
echo "Running bwa mem for ${SAMPLE_PREFIX} and piping to samtools for direct BAM conversion..."
# bwa mem outputs SAM to stdout, which is piped directly to samtools view
# samtools view -@ for threads (one less than total CPUs for BWA's main thread)
# -Sb - : Read SAM from stdin (-), output compressed BAM (-b), input is SAM (-S)
bwa mem -t $SLURM_CPUS_PER_TASK "$REFERENCE_GENOME" "$R1_FILE" "$R2_FILE" | \
  "$SAMTOOLS_EXE" view -@ $(($SLURM_CPUS_PER_TASK - 1)) -Sb - > "$OUTPUT_UNSORTED_BAM"
if [ $? -ne 0 ]; then
    echo "Error: BWA alignment or direct SAM to BAM conversion failed for ${SAMPLE_PREFIX}."
    exit 1
fi
echo "BWA alignment and BAM conversion for ${SAMPLE_PREFIX} complete."

# --- 2. Sort BAM ---
echo "Sorting BAM for ${SAMPLE_PREFIX}..."
# samtools sort -@ for threads, -o for output file
NUM_THREADS_SAMTOOLS=$(($SLURM_CPUS_PER_TASK - 1))
if [ "$NUM_THREADS_SAMTOOLS" -lt 1 ]; then NUM_THREADS_SAMTOOLS=1; fi # Ensure at least 1 thread

"$SAMTOOLS_EXE" sort -@ "$NUM_THREADS_SAMTOOLS" -o "$SORTED_BAM" "$OUTPUT_UNSORTED_BAM"
if [ $? -ne 0 ]; then
    echo "Error: Samtools sort for ${SAMPLE_PREFIX} failed."
    exit 1
fi
echo "Sorting for ${SAMPLE_PREFIX} complete."

# --- 3. Index Sorted BAM ---
echo "Indexing sorted BAM for ${SAMPLE_PREFIX}..."
"$SAMTOOLS_EXE" index "$SORTED_BAM"
if [ $? -ne 0 ]; then
    echo "Error: Samtools index for ${SAMPLE_PREFIX} failed."
    exit 1
fi
echo "Indexing for ${SAMPLE_PREFIX} complete."

# --- Optional: Clean up unsorted BAM ---
rm "$OUTPUT_UNSORTED_BAM"
echo "Removed unsorted BAM: $OUTPUT_UNSORTED_BAM"

echo "BWA Alignment and Samtools Processing Task for ${SAMPLE_PREFIX} Finished!"
echo "Date: $(date)"