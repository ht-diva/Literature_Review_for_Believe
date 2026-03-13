#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --job-name=decode2023_extract_alleles
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --array=0-99
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

source ~/.bashrc
conda activate ~/miniconda3/envs/gwasstudio

# File list
FILE_LIST="decode_2023_filenames_list.txt"
TOTAL_FILES=$(wc -l < "$FILE_LIST")

# Calculate files per job (ceiling division)
FILES_PER_JOB=$(( (TOTAL_FILES + 100 - 1) / 100 ))

# Determine start/end line for this array task
START=$(( SLURM_ARRAY_TASK_ID * FILES_PER_JOB + 1 ))
END=$(( START + FILES_PER_JOB - 1 ))

# Make sure END does not exceed total files
if [ "$END" -gt "$TOTAL_FILES" ]; then
    END=$TOTAL_FILES
fi

# Extract subset of files for this job
SUBSET=$(sed -n "${START},${END}p" "$FILE_LIST")

# Process each file sequentially
for FP in $SUBSET; do
    python process_single_file.py "$FP"
done

