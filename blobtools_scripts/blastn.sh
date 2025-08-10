#!/bin/bash -l
#SBATCH --job-name=BLASTn-array ##change to refelct sample running
#SBATCH --account=biocrusts
#SBATCH --partition=tier3
#SBATCH --time=0-20:00:00
#SBATCH --/shared/rc/microalgae/metagenomes/assemblies/MH_TOL5/blastn_logs/output=%x_%j.log
#SBATCH --/shared/rc/microalgae/metagenomes/assemblies/MH_TOL5/blastn_logs/error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4g
#SBATCH --mail-user=slack:    # Replace with your slack user
#SBATCH --mail-type=BEGIN,END,FAIL

# FIRST we need to split the assembly FASTA file and check the number of splits (Arrays):
## I use split_fast.py on my final.contigs.fa file then manually input that number below to set the correct number of arrays

## YOU ADD #SBATCH LINE BELOW TO MAKE A JOB ARRAY - (change 50 to number of split files)
#SBATCH --array=1-50

#Do not edit the echo sections
echo "All jobs in this array have:"
echo "- SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID}"
echo "- SLURM_ARRAY_TASK_COUNT=${SLURM_ARRAY_TASK_COUNT}"
echo "- SLURM_ARRAY_TASK_MIN=${SLURM_ARRAY_TASK_MIN}"
echo "- SLURM_ARRAY_TASK_MAX=${SLURM_ARRAY_TASK_MAX}"
echo "This job in the array has:"
echo "- SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"


# Load Genomics env - Blast+ module

spack env activate default-genomics-x86_64-25032001

# Run Blast comparisons 
##set directories as needed for split files and output

blastn -query /shared/rc/microalgae/metagenomes/assemblies/final.contigs${SLURM_ARRAY_TASK_ID}.fa -db /shared/rc/datasets/genomics/ncbi-blast-nt_2024-09-30/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads 15 -out /shared/rc/microalgae/metagenomes/assemblies/blastn_result/${SLURM_ARRAY_TASK_ID}-blastn.out
