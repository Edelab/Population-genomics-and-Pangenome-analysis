import numpy as np
import sys
from collections import defaultdict
from itertools import combinations

# Parsing user input
infile_list = sys.argv[1] #List of fasta
groups = sys.argv[2] # List of Group


#constant
OUTPUT_folder = "08_dnds"
INFILE_FOLDER = "07_Pnps/Cluster"
FASTA_FOLDER = "Clean_fasta/"
OUTPUT_FILE = "_DNDS.txt"

# Functions
def build_fasta_path(fasta_path, group = ""):
    return INFILE_FOLDER + group + "/"+ FASTA_FOLDER + fasta_path.strip()
    

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


def fast_dict_to_numpy(fasta_dict):
    seq_list = [fasta_dict[x] for x in fasta_dict]
    return np.array(seq_list)

def unequal_length(seq_len):
    return len(set(seq_len)) != 1

def null_length(seq_len):
    return 0 in set(seq_len)

def check_validity(codons_dict):
    seq_len = []
    for grp in codons_dict:
        seq_len.extend(np.apply_along_axis(len, arr= codons_dict[grp], axis = 1))


    if null_length(seq_len):
        return (False, "No sequence")

    if unequal_length(seq_len):
        return (False, "Unequal sequence length")

    return True, None


def parse_group(group_path):
    with open(group_path) as grp:
        return [x.strip() for x in grp.readlines()]


def parse_fasta(fast_path):
    seq_tab = defaultdict(list)
    
    file = open(fast_path, "r")
    indv = file.readline().strip()[1:]  
    seq = ""
    for line in file:
        if line.startswith(">"):
            seq_tab[indv].extend([seq.strip()[i:i+3] for i in range(0,len(seq),3)])
            indv = line.strip()[1:]
            seq = ""
            continue
        seq+= line.strip()

    seq_tab[indv].extend([seq.strip()[i:i+3] for i in range(0,len(seq),3)])
    return seq_tab


def fasta_reader(path, group1, group2):
   return {group1: fast_dict_to_numpy(parse_fasta(build_fasta_path(path, group1))),
              group2: fast_dict_to_numpy(parse_fasta(build_fasta_path(path, group2)))}


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

def merge_set(set1, set2):

    for i in [set1, set2]:
        if "NNN" in i:
            i.remove("NNN")
        if "---" in i:
            i.remove("---")
        if len(i) > 1:
            return {"NNN"}     
       
    return set1.union(set2)

def extract_mutations(codon_arr1, codon_arr2):
    set1 = np.apply_along_axis(set, arr= codon_arr1, axis = 0)
    set2 = np.apply_along_axis(set, arr= codon_arr2, axis = 0)

    mutations = [merge_set(codon, set2[x]) for x, codon in enumerate(set1)]
    return [x for x in mutations if len(x) >1]

def count_syn_nonsyn(codon_arr1, codon_arr2):
    mutation = extract_mutations(codon_arr1, codon_arr2)

    if len(mutation) == 0:
        return "0", "0"
    score = np.array(list(map(parse_codon_mut, mutation)))
    return str(np.sum(score[:,0])), str(np.sum(score[:,1])) 

def parse_sequence(path, group1, group2):

    print(f"parsing {path}")
    codon_dict = fasta_reader(path, group1, group2)
    usable, cause = check_validity(codon_dict)

    if not usable:
        print(f"skipping: {path.strip()}, when comparing: {group1}, {group2} reason: {cause}")
        return path.strip() + "\t-9\t-9"

    results = count_syn_nonsyn(codon_dict[group1],
                               codon_dict[group2])

    return path.strip() + "\t"+ "\t".join(results)



#parse sequence
list_file = open(infile_list, "r").readlines()
grouping = parse_group(groups)

group_pair = combinations(grouping, 2)

results = {(gr1,gr2): [parse_sequence(x, gr1, gr2) for x in list_file] for gr1, gr2 in group_pair}

for gr1, gr2 in results:
    with open(OUTPUT_folder + "/" + gr1 + "_" + gr2 + OUTPUT_FILE, "w") as outfile:
        outfile.write("Name\tDNS\tDS\n")
        outfile.write("\n".join(results[(gr1, gr2)]) + "\n")




