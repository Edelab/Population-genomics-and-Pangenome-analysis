# Long-read WGS preprocessing, assembly, annotation, and repeat masking pipeline for _P. brassicae_ isolates

## ONT read quality assessment with NanoPlot   

```bash

NanoPlot -t 32 --fastq sample_07_7.fastq.gz -o sample_7 --verbose -f png --tsv_stats --N50

```

## Reference genome mapping percentage with minimap2 and samtools

```bash

minimap2 -ax map-ont /directory/path/reference_genome_pb3A/Pb3A_genomic.fna sample_27_27.fastq.gz > sammple_27.sam

samtools view -bS sammple_27.sam | samtools sort -o sammple_27.sorted.bam

samtools index sammple_27.sorted.bam

samtools idxstats sammple_27.sorted.bam

```

## De-novo genome assembly using Hifiasm and Canu in ONT mode

## Hifiasm in ONT mode

### Download & Installation

```bash

# Download

wget https://github.com/chhylp123/hifiasm/releases/tag/0.25.0/hifiasm-0.25.0.tar.gz

# Installation

tar -xzvf hifiasm-0.25.0.tar.gz
cd hifiasm-0.25.0
make
./hifiasm --version

```

### Running Hifiasm

```bash

./hifiasm -t32 --ont -l0 --rl-cut --sc-cut --hg-size 25.5m --telo-m TTTTAGGG -o sample.asm --hom-cov auto sample.fastq.gz

## --rl-cut and --sc-cut, have been added to internally filter out low-quality reads before assembly.

# The FASTA file can be produced from GFA as follows:

awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa

```

### Running Canu

```bash
canu -p pb_canu -d canu_ont genomeSize=25.3m batMemory=60 -nanopore-raw raw_uncorrected_reads.fastq

#batMemory=60 to fix issue of canu using default batmemory=64MB if system has less bat memory

#Important code modifications to fix the issues consensus resuming and run CANU after correction to generate the final assembly.fa file  

useGrid=true \
# Give each consensus job 16 GB RAM and 8 threads
cnsMemory=16 \
cnsThreads=8 \
# Tell Canu to pass this to sbatch for each cns task
gridOptionsCns="--mem-per-cpu=2G -c 8 -p medium"

# The internal sbatch calls for each ctgcns/*.cns job will include only, and SLURM will accept them.

```

## Assembly contamination screening with NCBI Foreign Contamination Screen (FCS-GX)

```bash

singularity exec --bind /biodata:/app/db/gxdb/ --bind /directory/path/sample/fcs_clean:/sample-volume/ --bind /directory/path/sample/fcs_clean:/output-volume/ /directory/path/sample/fcs_clean/fcs-gx.sif python3 /app/bin/run_gx --fasta /sample-volume/qc56_assembly.fasta --out-dir /output-volume/ --gx-db /app/db/gxdb/gxdb --tax-id 37360

# cores = 30
# Memory = 512   

```

## Genome completeness assessment with BUSCO v5

```bash

# sed command to remove the * from the end of the FASTA file
sed -E '/^>/!s/\*$//' input.fasta > output.fasta

# Load modules
module load miniconda/3-py3.10 busco/5.8.2

# Initialize conda environment

eval "$(conda shell.bash hook)"
conda activate busco-5.8.2

busco -i /path/to/genome_assembly/contigs.fasta -o output -l /path/to/lineage/eukaryota_odb10 –-mode genome -c 6

# Use **--force** if need to overwrite the existing directory)

# If the lineage needs to be downloaded, simply write the **-l eukaryota_odb10**; this will download the lineage and can be used for further runs

```

## Genome completeness assessment with Compleasm

```bash

python -m compleasm run -a sample_hifiasm_asm.fa -o sample_raw_assembly -t 20 -l eukaryota

```

## Hapo_G polishing with Illumina short reads

```bash

module load bwa/0.7.13 htslib/1.8 samtools/1.13 python/3.7 hapog/1.2
/prg/hapog/1.2/hapog.py --genome sample_assembly_hifiasm.fa --pe1 sample_forward_paired.fastq.gz --pe2 sample_reverse_paired.fastq.gz -o hapoG -t 10 -u

```

## Repeat element prediction and masking with RepeatModeler2 and RepeatMasker of the polished genome assemblies

```bash

# RepeatModeler
module load RepeatModeler/2.0.1
BuildDatabase -name qc56_genome_db -engine ncbi sample_canu_hapoG.fa

RepeatModeler -database sample_genome_db -pa 8 -LTRStruct > repeatmodeler.log

# Repeat masking

RepeatMasker/4.1.2
RepeatMasker -pa 8 -lib /directory/path/consensi.fa.classified -gff -xsmall -dir repeatmasker_output sample/polished_round_2.fa

```

## Transposable element annotation with EarlGrey (v6)

## Installations

```bash

# With conda
conda create -n earlgrey -c conda-forge -c bioconda earlgrey=6.0.3     # Recent version

conda activate earlgrey
earlGrey --help

```
## Running TE annotations

```bash

earlGrey -g genome_assembly.fa -t 25 -c yes -d yes -r 2759 -s plasmodiophoraBrassicae -o ./earlGreyOutputs

```

## Gene prediction with BRAKER3 using protein and RNA-seq evidence

```bash

#!/bin/bash

#SBATCH -J annotation
#SBATCH -o Pb_anno.out
#SBATCH -e Pb_anno.err
#SBATCH -c 48
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@ulaval.ca
#SBATCH --time=1-00:00
#SBATCH --mem=100G

# Load necessary modules
#module load singularity

# Define paths
export AUGUSTUS_CONFIG_PATH=/directory/path/config
WORKDIR=/directory/path/
GENOME=/directory/path/polished_assembly.fa
PROTEINS=proteins.fa
RNASEQ_DIR=/directory/path/rnaseq

# Ensure working directory exists
mkdir -p $WORKDIR

# Run BRAKER3
singularity exec -B ${PWD}:${PWD} -B /directory/path/config:/opt/Augustus/config braker3.sif braker.pl \
    --AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH \
    --useexisting \
    --genome=$GENOME \
    --species=isolate_name \
    --prot_seq=$PROTEINS \
    --rnaseq_sets_ids=S001A93_ID1,S001E7A_ID2,21dpi_ID3,7dpi_ID4 \
    --rnaseq_sets_dirs=$RNASEQ_DIR \
    --workingdir=$WORKDIR \
    --threads=48 \
    --busco_lineage=eukaryota_odb10

```

## Extract longest isoform per gene from BRAKER GTF based on CDS coordinate length 

```python

#!/usr/bin/env python3

# ============================================================================================
# comparing_fa_gtf_longest_isoforms.py
# Keeps only the longest isoform per gene from BRAKER3 output GTF based on CDS coordinate length
# python comparing_fa_gtf_longest_isoforms.py braker.gtf > braker_longest_isoforms.gtf
# ============================================================================================


import csv
import sys
import re
import argparse


def extractFeature(text, feature):
    regex = feature + ' "([^"]+)"'
    result = re.search(regex, text)
    if result:
        return result.groups()[0]
    else:
        return None


def computeLengths(input):
    transcriptLengths = dict()
    for row in csv.reader(open(input), delimiter='\t'):
        if len(row) == 0 or row[0].startswith('#'): continue
        if (row[2] == 'CDS'):
            gene = extractFeature(row[8], 'gene_id')
            transcript = extractFeature(row[8], 'transcript_id')
            if not gene or not transcript:
                continue
            if gene not in transcriptLengths:
                transcriptLengths[gene] = dict()
            if transcript not in transcriptLengths[gene]:
                transcriptLengths[gene][transcript] = 0
            transcriptLengths[gene][transcript] += int(row[4]) - int(row[3])
    return transcriptLengths


def getLongestTranscript(transcriptLengths):
    longestTranscripts = dict()
    for gene in transcriptLengths:
        max = 0
        longestTranscript = ""
        for transcript in transcriptLengths[gene]:
            length = transcriptLengths[gene][transcript]
            if (length > max):
                max = length
                longestTranscript = transcript
        longestTranscripts[gene] = longestTranscript
    return longestTranscripts


def printLongest(input, longestTranscripts):
    for row in csv.reader(open(input), delimiter='\t'):
        if len(row) == 0 or row[0].startswith('#'): continue
        gene = extractFeature(row[8], 'gene_id')
        transcript = extractFeature(row[8], 'transcript_id')
        if not gene or not transcript:
            continue
        if (longestTranscripts[gene] == transcript):
            print('\t'.join(row))


def main():
    args = parseCmd()
    transcriptLengths = computeLengths(args.input)
    longestTranscripts = getLongestTranscript(transcriptLengths)
    printLongest(sys.argv[1], longestTranscripts)


def parseCmd():

    parser = argparse.ArgumentParser(description='Print longest isoforms in a \
                                     gtf file')

    parser.add_argument('input', type=str,
                        help='Input gtf file')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()

```
## Attribution
Script `comparing_fa_gtf_longest_isoforms.py` adapted from Tomas Bruna,
Georgia Institute of Technology (2019).
Original: https://github.com/gatech-genemark/ProtHint/blob/master/bin/print_longest_isoform.py

