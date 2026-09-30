#!/bin/bash

#SBATCH -D /path/to/directory/array_job
#SBATCH -J repeat_annotation
#SBATCH -o repeat_logs/repeat_%A_%a.out
#SBATCH -e repeat_logs/repeat_%A_%a.err
#SBATCH -c 8
#SBATCH -p small
#SBATCH --time=2-00:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@ulaval.ca
#SBATCH --array=0-13  # Adjust according to number of samples

# Load required modules
module load RepeatModeler/2.0.1
module load RepeatMasker/4.1.2

# Read sample from file
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

# Set paths
workdir="/path/to/directory/array_job/${sample}"
genome="${workdir}/hapog_3/hapog_results/${sample}_hapog_3.fasta"

# Set output directories
rm_outdir="${workdir}/repeatmodeler_output"
mask_outdir="${workdir}/repeatmasker_output"

# Create output dirs
mkdir -p "$rm_outdir" "$mask_outdir"
cd "$rm_outdir"

# Run RepeatModeler
echo "Running RepeatModeler for $sample ..."
BuildDatabase -name ${sample}_genome_db -engine ncbi "$genome"

RepeatModeler -database ${sample}_genome_db \
              -pa 8 \
              -LTRStruct \
              > repeatmodeler.log

# Find classified library output
libfile=$(find "$rm_outdir" -name "consensi.fa.classified" | head -n 1)

# Run RepeatMasker if library was created
if [[ -f "$libfile" ]]; then
    echo "Running RepeatMasker for $sample using library: $libfile"
    RepeatMasker -pa 8 \
                 -lib "$libfile" \
                 -gff \
                 -xsmall \
                 -dir "$mask_outdir" \
                 "$genome"
else
    echo "ERROR: RepeatModeler library not found for $sample"
    exit 1
fi

echo "RepeatModeler and RepeatMasker finished for $sample"
