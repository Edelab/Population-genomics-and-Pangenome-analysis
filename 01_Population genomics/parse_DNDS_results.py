#!/bin/python3

"""
Script puting together result from
1)sumstat_seq
2_DNDS

"""
# imports
import numpy as np
from collections import defaultdict
import sys 
from itertools import combinations


# Constants
CLUSTER_PREFIX = "Cluster"
PNPS_FOLDER = "07_Pnps/"
DNDS_FOLDER = "08_dnds/"
DNDS_SUFFIX = "_DNDS.txt"
SUMSTATS_SUFFIX = ".CDS.sumstats"

groups = sys.argv[1]
output_file = sys.argv[2] #"POP_withoutquantiles_completesumstats"

# Functions
def parse_group(group_path):
    with open(group_path) as grp:
        return [x.strip() for x in grp.readlines()]

def build_sumstats_path(group):
    return PNPS_FOLDER + CLUSTER_PREFIX + group + "/" + CLUSTER_PREFIX + group + SUMSTATS_SUFFIX 

def build_dnds_path(gr1, gr2):
    return DNDS_FOLDER + gr1 + "_" + gr2 + DNDS_SUFFIX



def parsing_input(path: str):
    parsed_file = defaultdict(lambda :defaultdict(float))
    infile = open(path, "r")
    l = infile.readline().strip().split("\t")
    header = {x.lower(): y for y, x in enumerate(l)} 
    for line in infile:
        l = line.strip().split("\t")
        name = l[header["name"]].split(".")[0]
        metrics = {x: float(l[header[x]]) for x in header if x != "name"}
        parsed_file[name] = metrics
    return parsed_file



### MODIFY THIS TO CALCULATE ESTIMATES
def merging_dataset(sumstats1, sumstats2, dnds) -> dict:
    merged_data = defaultdict(list)
    for gene in dnds:

        n = (sumstats1[gene]["n"] + sumstats2[gene]["n"])/2
        size = (sumstats1[gene]["size"] + sumstats2[gene]["size"])/2

        if size == 0 or dnds[gene]["dns"] == -9 :
            merged_data[gene] = {"n" : n, "size" : size, "DN": -999, "DS" :  -999 }
            continue
        
        ss = (sumstats1[gene]["nss"] + sumstats2[gene]["nss"])/2
        Nss = size - ss

        merged_data[gene] = { "n" : n, "size" : size, "DN": dnds[gene]["dns"] / Nss, "DS" :  dnds[gene]["ds"] / ss }

        
    return merged_data


def calculte_dnds_ratio(data):
    return None

def simplify_merge(data):
    return None
    

def conv_to_np(data):
    header = list(data[next(iter(data))].keys())
    data_array = np.zeros((len(data), int(len(header))))
    gene_order = [0]* len(data)

    for x, gene in enumerate(data):
        gene_order[x] = gene
        data_array[x,:len(header)] = [data[gene][y] for y in header]
    return {"gene": gene_order, "header": header, "metrics": data_array.astype(float)}


  
 

grouping = parse_group(groups)

group_pair = list(combinations(grouping, 2))

sumstats = {gr:  parsing_input(build_sumstats_path(gr)) for gr in grouping }
dnds = {(gr1, gr2):  parsing_input(build_dnds_path(gr1, gr2)) for gr1, gr2 in group_pair }

merged_dataset = {}
for gr1, gr2 in group_pair:
    tmp = merging_dataset(sumstats[gr1], sumstats[gr2], dnds[(gr1, gr2)])
    merged_dataset[(gr1, gr2)] = conv_to_np(tmp)


##not working now
##outputing files
for gr1, gr2 in group_pair:
    out = open(output_file + f"{gr1}_{gr2}.DNDNS.final", "w")
    out.write("name" + "\t" + "\t".join(merged_dataset[(gr1,gr2)]["header"]) + "\n")
    final_str = list(zip(merged_dataset[(gr1, gr2)]["gene"], ["\t".join(list(map(str,x))) for x in merged_dataset[(gr1, gr2)]["metrics"]]))
    out.write("\n".join(["\t".join(x) for x in final_str]))




