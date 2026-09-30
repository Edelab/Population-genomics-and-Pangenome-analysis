
## Phylogenetic tree construction from single-copy orthologue alignment using IQ-TREE2

```bash

#!/bin/bash

outgroup="02_outgroup_genes_renamed"
sample="03_samples_fasta"
concat="04_concatenated_fasta"
prot="05_protein_fasta"
aligned="06_aligned_protein"
aligned_filtered="07_cleaned_alignement"
reversed="08_reversed_alignement"
# Merge outgroup and sample in one fasta for each gene
for i in $(ls $outgroup/)
do
        gene=${i/.fasta}
        echo $gene
        cat $outgroup/$gene.fasta $sample/$gene.fst >$concat/$i
        python 01_scripts/00_utils/Fasta_to_prot.py $concat/$i $prot/${i/.fasta}
        muscle -align $prot/${i/.fasta}.prot -output $aligned/${i/.fasta}.prot
done

# Filtering alignement
Rscript 00_utils/selection_gene.R
python3 00_utils/filter_alignement.py $aligned blacklisted.gene $aligned_filtered

for fasta in $aligned_filtered/*.prot
    do
    python3 00_utils/alignement_transferer.py $fasta $concat $reversed
    done


# Concatenate genes sequence based on sample name
seqkit concat $reversed/* -o 09_concatenate_sequences/concatenated.fa -j 10

# Filter non polymorphic positions
python3 00_utils/filter_fasta.py 09_concatenate_sequences/concatenated.fa 09_concatenate_sequences/polymorphic.concatenated.fasta Spongospora

# Remove NA
python3 00_utils/filter_merged_alignment.py 09_concatenate_sequences/polymorphic.concatenated.fasta 09_concatenate_sequences/no_missing_polymorphic.concatenated.fasta 0

# run iqtree

iqtree2 -s 09_concatenate_sequences/no_missing_polymorphic.concatenated.fasta -m MFP+ASC -B 1000

```

1. Nucleotide to protein translation of concatenated FASTA sequences for phylogenetic alignment

```python

import sys
filename = sys.argv[1]
outprefix=sys.argv[2]

def translate_dna(sequence):

        codontable = {
                    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
                    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
                    '---': '-'}

        proteinsequence = ''
        sequencestart = sequence[0:]
        cds = str(sequencestart[:len(sequence)+3])
        for n in range(0,len(cds),3):
                if cds[n:n+3].upper() in codontable:
                        proteinsequence += codontable[cds[n:n+3].upper()]
                else:
                        proteinsequence += "X"

        return proteinsequence

                                                                                
infilename=open(filename)
outfilenameprot=outprefix + ".prot"
outfilenameprotopen=open(outfilenameprot,"w")


for line in infilename:
    if line[0] == ">":
        outfilenameprotopen.write(line)
    else:
        line2print=translate_dna(line.strip())+"\n"
        outfilenameprotopen.write(line2print)

infilename.close()

```

2. Gene removal under selection based on dN/dS and pN/pS outlier thresholds across three clusters

```R

###Gene filtering for selection

dnds_1_2 <- read.table("dnds/1_2.DNDNS.final", h= T)
dnds_1_2$ratio_1_2 <- dnds_1_2$DN/dnds_1_2$DS
dnds_1_2 <- dnds_1_2[dnds_1_2$DS >0 & dnds_1_2$n >20,-c(2,4,5)]

dnds_2_3 <- read.table("dnds/2_3.DNDNS.final", h= T)
dnds_2_3$ratio_2_3 <- dnds_2_3$DN/dnds_2_3$DS
dnds_2_3 <- dnds_2_3[dnds_2_3$DS >0 & dnds_2_3$n>20 ,-c(2,4,5)]


dnds_1_3 <- read.table("dnds/1_3.DNDNS.final", h= T)
dnds_1_3$ratio_1_3 <- dnds_1_3$DN/dnds_1_3$DS
dnds_1_3 <- dnds_1_3[dnds_1_3$DS >0 & dnds_1_3$n > 20, -c(2,4,5)]


DNDS <- merge(dnds_1_2, dnds_1_3, by = c("name", "size"), sort = FALSE, all = TRUE)
DNDS <- merge(DNDS, dnds_2_3, by = c("name", "size"), sort = FALSE, all = TRUE)

limit_dnds <- quantile(c(DNDS$ratio_1_2, DNDS$ratio_1_3, DNDS$ratio_2_3), 0.95, na.rm = T)


bad_genes_dnds <- c()

for(i in 1:nrow(DNDS)){
  a <- DNDS$ratio_1_2[i] > limit_dnds
  b <- DNDS$ratio_1_3[i] > limit_dnds
  c <- DNDS$ratio_2_3[i] > limit_dnds
  
  if( any(a,b,c, na.rm=T)){
    bad_genes_dnds <- c(bad_genes_dnds, DNDS$name[i])
  }
}


PNPS_1 <- read.table("pnps/Cluster1/Cluster1.final.cleaned", h=T)
PNPS_1[is.na(PNPS_1$ps),] <- -9
PNPS_1 <- PNPS_1[PNPS_1$ps >0 & PNPS_1$size >0,]
PNPS_1$PNPS <- PNPS_1$pns/PNPS_1$ps

PNPS_2 <- read.table("pnps/Cluster2.final.cleaned", h=T)
PNPS_2[is.na(PNPS_2$ps),] <- -9

PNPS_2 <- PNPS_2[PNPS_2$ps >0 ,]
PNPS_2$PNPS <- PNPS_2$pns/PNPS_2$ps


PNPS_3 <- read.table("pnps/Cluster3.final.cleaned", h=T)
PNPS_3[is.na(PNPS_3$ps),] <- -9

PNPS_3 <- PNPS_3[PNPS_3$ps >0,]
PNPS_3$PNPS <- PNPS_3$pns/PNPS_3$ps

PNPS <- merge(PNPS_1, PNPS_2, by = c("name", "size"), sort = FALSE, all = TRUE)
PNPS <- merge(PNPS, PNPS_3, by = c("name", "size"), sort = FALSE, all = TRUE)


limit_pnps <- quantile(c(PNPS_1$PNPS, PNPS_2$PNPS, PNPS_3$PNPS), 0.95, na.rm = T)
bad_genes_pnps <- c()

for(i in 1:nrow(PNPS)){
  a <- PNPS$PNPS[i] > limit_pnps
  b <- PNPS$PNPS.x[i] > limit_pnps
  c <- PNPS$PNPS.y[i] > limit_pnps
  
  if( any(a,b,c, na.rm=T)){
    bad_genes_pnps <- c(bad_genes_pnps, PNPS$name[i])
  }
}


bad_genes <- unique(c(bad_genes_dnds, bad_genes_pnps))
write.table(x = bad_genes,
            file="blacklisted.gene",quote =  FALSE,row.names = FALSE, col.names = FALSE)
```


3. Filter protein alignments by checking start codon stop codon and gap threshold

```python

#!/bin/env python


## Imports
import sys
import glob
import shutil

#Functions
def parse_blacklist(path):
    return list(map(lambda x: x.strip(), open(path, "r").readlines()))

def extract_gene_name_from_path(path):
    return path.split("/")[-1].split(".")[0]

def Check_blacklist(path, blacklist):
    return extract_gene_name_from_path(path) in blacklist

def parse_ali(path):
    alif = open(path, "r")
    ali = {}
    for line in alif:
        if line.startswith(">"):
            SN = line.strip()[1:]
            ali[SN] = ""
            continue
        
        ali[SN] += line.strip()
    alif.close()    
    return ali

def check_M(ali):
    return not any([ali[x][0] not in ["M", "X"] for x in ali])


def check_stop(ali):
    return any(["_" in ali[x][:-1] for x in ali])

def get_gap(ali):
    return max([len(ali[x]) - len(ali[x].replace("-", "")) for x in ali])

def build_output_path(file, outputfolder):
    return outputfolder + "/" + extract_gene_name_from_path(file) + ".prot"

def write_ali(ali, path):
    ali_line = ""
    for ind in ali:
        ali_line += f">{ind}\n{ali[ind]}\n"

    file = open(path, "w")
    file.write(ali_line)
    file.close()

#Main
def main():
    ali_folder, blacklist_path, output_folder = sys.argv[1:4]
    ali_files = glob.glob(ali_folder + "/*.prot")
    blacklist = parse_blacklist(blacklist_path)
    for file in ali_files:
        if Check_blacklist(file, blacklist):
            print(f"Skipping {file}: Blacklisted")
            continue
        ali = parse_ali(file)
        if not check_M(ali):
            print(f"Skipping {file}: No start codon")
            continue
        if check_stop(ali):
            print(f"Skipping {file}: intermediate stop codon")
            continue
        if get_gap(ali) > 100:
            print(f"Skipping {file}: Too many gaps")
            continue
        
        write_ali(ali, build_output_path(file, output_folder))



if __name__ == "__main__":
    main()

```

4. Transfer protein alignment gaps to nucleotide sequences for codon-aware alignment

```python

#!/bin/env python

import sys


#user variables

prot_al_path = sys.argv[1]
nt_al_path = sys.argv[2]
out = sys.argv[3]

# function

def parse_seq(path):
    seqs = dict()
    with open(path, "r") as seq:
        ind = ""
        for line in seq:
            l = line.strip()
            if l.startswith(">"):
                ind = l[1:]
                seqs[ind] = ""
                continue
            seqs[ind] += l
    return seqs

def write_fasta(path, seqdict):
    with open(path, "w") as out:
        for ind in seqdict:
            out.write(f">{ind}\n{seqdict[ind]}\n")




# parsing seqs

prot_seq = parse_seq(prot_al_path)
nuc_seq = parse_seq(nt_al_path)


# parse sequence_shift:
seq_shift = dict()

for ind in prot_seq:
    seq_shift[ind] = []
    count = False
    counter = 0
    start = -9

    for x, aa in enumerate(prot_seq[ind]):
        if aa == "-":
            if not count:
                count = True
                start = x
            counter += 1

        if aa != "-":
            if count:
                count = False
                seq_shift[ind].append((start, counter))
                counter = 0
                start = -9
    
    if count:
        seq_shift[ind].append((start, counter))


# Reconstruct NT alignement from sequence shift
newseq = dict()
for ind in seq_shift:
    tmpseq = nuc_seq[ind]
    for pos, count in seq_shift[ind]:
        tmpseq = tmpseq[:(pos * 3 )] + "-" * 3 * count + tmpseq[(pos*3):]
    newseq[ind] = tmpseq

write_fasta(out, newseq)

```

5. Filter concatenated FASTA to polymorphic sites excluding N and gap characters

```python

#!/bin/python3

from collections import defaultdict
import sys


fasta= sys.argv[1]
output=sys.argv[2]
if len(sys.argv) == 4:
    outgroup = sys.argv[3]
else:
    outgroup = ""
outgroup_seq = ""


def is_polymorphic(pos):
    seq = [fastadict[x][pos] for x in fastadict]
    if pos == 4116182:
        print(seq)

        print(set(seq))
    if len(set([x for x in seq if x not in ["N", "-"]])) <= 1:
        return None
    return pos

print("Parsing sequences...")
with open(fasta, "r") as infile:
    fastadict = defaultdict(str)

    for line in infile:
        l = line.strip()
        if line.startswith(">"):
            ind = l[1:]
            continue
        if ind == outgroup:
            outgroup_seq +=l
            continue
        fastadict[ind] += l


print("Finding polymoprhic sites")


pos = [x for x in map(is_polymorphic, range(len(outgroup_seq))) if x is not None]
print(len(pos))

print("Writing results")
with open(output, "w") as out:
    for i in fastadict:

        out.write(">" + i + "\n")
        out.write("".join([fastadict[i][x] for x in pos]) + "\n")
    
    out.write(">" + outgroup + "\n")
    out.write("".join([outgroup_seq[x] for x in pos]) + "\n")
```

6. Filter merged alignment sequences exceeding missing data threshold before phylogenetic analysis

```python

#!/bin/env python3

import sys

infile, outfile, perc_missing = sys.argv[1:4]

def parse_ali(path):
    alif = open(path, "r")
    ali = {}
    for line in alif:
        if line.startswith(">"):
            SN = line.strip()[1:]
            ali[SN] = ""
            continue

        ali[SN] += line.strip()
    alif.close()
    return ali


ali = parse_ali(infile)
perc_missing = float(perc_missing)

with open(outfile, "w") as out:
    for ind in ali:
        if ali[ind].count("N") / len(ali[ind]) >= perc_missing:
            print(f"filtering {ind} out... {ali[ind].count("N") / len(ali[ind])}% missing data")
            continue
        out.write(f">{ind}\n{ali[ind]}\n")

```
        

