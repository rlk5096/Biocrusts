#!/bin/bash -l
#SBATCH --job-name=undef_blast-array        ##change to reflect job 
#SBATCH --account=biocrusts
#SBATCH --partition=tier3
#SBATCH --time=0-20:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4g
#SBATCH --mail-user=slack:      # Add your slack user
#SBATCH --mail-type=BEGIN,END,FAIL

##after extracting no-hit fasta, parse and manually add resulting number of splits to set array number 
# split_fasta.py

## YOU ADD #SBATCH LINE BELOW TO MAKE A JOB ARRAY (chagne 50 to number of split fasta files)
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

##running diamond blastx - this is from a version installed in the home directory
##change query and out paths - be careful in making sure directories specifed are already created or add makedir command to script
/home/rlk5096/diamond blastx \
    --db /shared/rc/datasets/genomics/ncbi-blast-nr_2024-10-01/nr \
    --query /shared/rc/microalgae/metagenomes/assemblies/no_hitsta${SLURM_ARRAY_TASK_ID}.fa \
    --outfmt 6 \
    --max-target-seqs 1 \
    --evalue 1e-25 \
    --threads 11 \
    --out /shared/rc/microalgae/metagenomes/assemblies/${SLURM_ARRAY_TASK_ID}_nohit_match.tsv
