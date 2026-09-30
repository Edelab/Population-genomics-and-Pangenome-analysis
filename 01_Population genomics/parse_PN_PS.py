import numpy as np
import sys
from collections import defaultdict


# Parsing user input
infile_list = sys.argv[1]
prefix = sys.argv[2]

#constant
INFILE_FOLDER = "Clean_fasta"
OUTPUT_FILE = "_PNPS.txt"
# Functions

def translate_condon(codon):
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
                    }
        return codontable[codon] if codon in codontable else "X"


def parse_fasta(fast_path):
    seq_tab = list()
    file = open(fast_path, "r")
    file.readline()
    seq = ""
    for line in file:
        if line.startswith(">"):
            seq_tab.append([seq[i:i+3] for i in range(0,len(seq),3)])
            seq = ""
            continue
        seq += line.strip()
    seq_tab.append([seq[i:i+3] for i in range(0,len(seq),3)])
    return np.array(seq_tab)


def extract_mutations(codon_arr):
    codon_set = np.apply_along_axis(set, arr= codon_arr, axis = 0)
    [x.discard("NNN") for x in codon_set]
    return [x for x in codon_set if len(x) >1]

def parse_codon_mut(codon):
    codon = np.array([list(x) for x in codon])
    mutation = dict()
    NS, S = defaultdict(int), defaultdict(int)

    for i in range(3):
         mutation[i] = set(codon[:,i])
    SNP = {x: mutation[x] for x in mutation if len(mutation[x]) >1}
    
    possible_codon = []
    for one in mutation[0]:
        for two in mutation[1]:
            for three in mutation[2]:
                possible_codon.append([one, two, three])
    
    for ref in possible_codon:
        ref_codon = translate_condon("".join(ref))
        NS["".join(ref)] = 0
        S["".join(ref)] = 0

        for i in SNP:
           muted = [np.concatenate((ref[:i], [x], ref[i+1:]), axis = 0) for x in mutation[i]]
           muted_codon = [translate_condon("".join(x)) for x in muted]
           if any(x != ref_codon for x in muted_codon):
               NS["".join(ref)] += 1
           else:
               S["".join(ref)] += 1
    return (min(NS.values()), len(SNP) - min(NS.values()))


def count_syn_nonsyn(codon_arr):
    mutation = extract_mutations(codon_arr)
    if len(mutation) == 0:
        return "0", "0"

    score = np.array(list(map(parse_codon_mut, mutation)))
    return str(np.sum(score[:,0])), str(np.sum(score[:,1])) 

def parse_sequence(path):
    codon_array = parse_fasta(INFILE_FOLDER + "/" + path.strip())
    if codon_array.shape[1] == 0:
        return path.strip() + "\t0\t0"
    return path.strip() + "\t"+ "\t".join(count_syn_nonsyn(codon_array))



#parse sequence
list_file = open(infile_list, "r").readlines()

result = map(parse_sequence, list_file)
with open(prefix + OUTPUT_FILE, "w") as outfile:
    outfile.write("Name\tPNS\tPS\n")
    outfile.write("\n".join(result))




