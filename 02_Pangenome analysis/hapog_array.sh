#!/bin/bash

#SBATCH -D /path/to/directory/array_job
#SBATCH -J hapog_array
#SBATCH -o logs/hapog_%A_%a.out
#SBATCH -e logs/hapog_%A_%a.err
#SBATCH -c 6
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@ulaval.ca
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G
#SBATCH --array=0-13  # Adjust based on the number of samples in samples.txt

# Load modules
module load bwa/0.7.13 htslib/1.8 samtools/1.13 python/3.7 hapog/1.2

# Read the sample name from file using SLURM_ARRAY_TASK_ID
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" fixed_samples.txt)

# Normalize sample name for case-sensitivity if needed
lower=$(echo "$sample" | tr '[:upper:]' '[:lower:]')

# Set paths
basedir="${sample}"
genome="${basedir}/${sample}_trimmed_hifiasm.fasta"
pe1="${basedir}/${sample}_forward_paired.fastq.gz"
pe2="${basedir}/${sample}_reverse_paired.fastq.gz"
output_dir="${basedir}/hapog_1"

# Run Hapo-G
echo "Running Hapo-G for $sample ..."
/prg/hapog/1.2/hapog.py \
    --genome "$genome" \
    --pe1 "$pe1" \
    --pe2 "$pe2" \
    -o "$output_dir" \
    -t 6 \
    -u
