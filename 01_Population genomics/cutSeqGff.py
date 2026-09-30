#! usr/bin/env python3

from __future__ import print_function
import sys
from collections import defaultdict


# cut a sequences using a GFF

def read_fasta(fasta_path):
    sample_id, sequence = None, []
    for l in fasta_path:
        l = l.rstrip()
        if l.startswith(">"):
            if sample_id:
                yield sample_id, ''.join(sequence)
            sample_id, sequence = l, []
        else:
            sequence.append(l)
    if sample_id:
        yield sample_id, ''.join(sequence)
    
 
def revcomp(sequence):
    return sequence.translate(sequence.maketrans('ACGTacgtRYMKrymkVBHDvbhdNn-', 'TGCAtgcaYRKMyrkmBVDHbvdhNn-'))[::-1]


def find_atr(attributes_str, atr):
    gene_name = [x.split("=")[1] for x in attributes_str.split(";") if x.startswith(atr)]
    return gene_name[0] if len(gene_name) >0 else ""


def sum_of_diff(exon_list):
    return abs(sum(exon_list[0]) - sum(exon_list[1]))


def extract_longest(transcript):
        size =  {tr :sum_of_diff([transcript[tr]["beg"], transcript[tr]["end"]]) for tr in transcript}
        
        maxsize = 0
        longest_tr = ""
        for tr in size:
            if size[tr] > maxsize:
                longest_tr = tr
                maxsize = size[tr]
        return longest_tr


def parse_cds(ref, transcript ):
    longest_tr = transcript[extract_longest(transcript)]
    output_str = ""
    
    coord = list(zip(longest_tr["beg"], longest_tr["end"]))
    coord.sort()

    if longest_tr["revc"][0] == "+":
        parser_seq = lambda x: x

    else:
        parser_seq = revcomp

    for name in ref:
        output_str += f">{name}\n"
        tmp = ""
        for b in coord:
            tmp += ref[name][(b[0] - 1):b[1]]
        output_str += parser_seq(tmp) + "\n"
    return {"transcript": longest_tr["geneName"], "seq": output_str}


# Main

fp = open(sys.argv[1])
GFFfile = open(sys.argv[2]) 
IdSeqIn = str(sys.argv[3]) 
output = str(sys.argv[4])

CDS_pos = defaultdict(lambda: defaultdict(list))
gene_tr = defaultdict(list)
ref_sequences = dict()



for name, seq in read_fasta(fp):
    ref_sequences[name.lstrip('>')] = seq


# Parsing CDS's
for x,line in enumerate(GFFfile):
    line = line.rstrip()
    if line[0] == "#":
        continue

    arrline = line.split()
    if arrline[0] != IdSeqIn:
        continue
    
    if arrline[2] != "CDS":
        continue
    
    attributes = arrline[8]
    geneName = find_atr(attributes, "gene_id")
    
    if geneName == "": #Some special CDS don't have name such as V_gene_segments. Could if using parent name instead?
        print(f" no name in {attributes}")
        continue
    
    tr = find_atr(attributes, "Parent")
    if tr == "":
        print(f"no transcript name in {attributes}")
        continue

    if tr not in gene_tr[geneName]:
        gene_tr[geneName].append(tr)
    
    CDS_pos[tr]["geneName"] = geneName
    CDS_pos[tr]["revc"].append(arrline[6])
    CDS_pos[tr]["beg"].append(int(arrline[3]))
    CDS_pos[tr]["end"].append(int(arrline[4]))
    CDS_pos[tr]["phase"].append(int(arrline[7]))


for gene in gene_tr:
    transcripts =  {tr: CDS_pos[tr] for tr in gene_tr[gene]}
    longest_tr = parse_cds(ref_sequences, transcripts)
    fasta = open(output + "/" + longest_tr["transcript"] + ".fst", "w")
    fasta.write(longest_tr["seq"])

