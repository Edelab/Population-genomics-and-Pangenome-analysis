# Custom python script to compare the .fa and .gtf longest isoforms to avoid discrapencies

from Bio import SeqIO
import pandas as pd
import re

# Paths to your files
FASTA_FILE = "6c_proteins_longest_isoform.fa"
GTF_FILE = "6c_braker_longest_isoform.gtf"

def get_fasta_ids(fasta_file):
    ids = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        ids.add(record.id)
    return ids

def get_gtf_ids(gtf_file):
    ids = set()
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            attr = cols[8]
            match = re.search(r'transcript_id "?([^\s";]+)"?', attr)
            if match:
                ids.add(match.group(1))
    return ids

# Extract ids
fasta_ids = get_fasta_ids(FASTA_FILE)
gtf_ids = get_gtf_ids(GTF_FILE)

# Comparisons
only_in_fasta = sorted(fasta_ids - gtf_ids)
only_in_gtf = sorted(gtf_ids - fasta_ids)
in_both = fasta_ids & gtf_ids

# Summary
print(f"Total IDs in FASTA: {len(fasta_ids)}")
print(f"Total IDs in GTF: {len(gtf_ids)}")
print(f"IDs present in both: {len(in_both)}")

if not only_in_fasta and not only_in_gtf:
    print("\n✅ All isoforms match between filtered FASTA and filtered GTF.")
else:
    print("\n❗ There are discrepancies. See below:")

# Show discrepancies as dataframes for easy review
if only_in_fasta:
    print("\nIsoforms only in FASTA but not in GTF:")
    display(pd.DataFrame(only_in_fasta, columns=["Isoform_ID"]))
if only_in_gtf:
    print("\nIsoforms only in GTF but not in FASTA:")
    display(pd.DataFrame(only_in_gtf, columns=["Isoform_ID"]))
