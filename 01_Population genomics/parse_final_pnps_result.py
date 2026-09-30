#!/bin/python3

"""
Script puting together result from
1)sumstat_seq
2)fourfold script
3)parse_PnPs scrit
and estimating Tajima's D
"""
# imports
import numpy as np
from collections import defaultdict
import sys 

# Constant path
sumstat_file = sys.argv[1]
fourfold_file = sys.argv[2] # "POP_withoutquantiles.4fold.CDS.sumstats"
PnPs_file = sys.argv[3] # "POP_withoutquantiles.CDS.sumstats_PNPS"
output_file = sys.argv[4] #"POP_withoutquantiles_completesumstats"

# Functions
def parsing_input(path: str):
    parsed_file = defaultdict(lambda :defaultdict(float))
    infile = open(path, "r")
    l = infile.readline().strip().split("\t")
    header = {x.lower(): y for y, x in enumerate(l)} 
    for line in infile:
        l = line.strip().split("\t")
        name = l[header["name"]].split(".")[0]
        metrics = {x: l[header[x]] for x in header if x != "name"}
        parsed_file[name] = metrics
    return parsed_file

    
def merging_dataset(stat: dict, fold4: dict, pnps: dict) -> dict:
    merged_data = defaultdict(list)
    for gene in stat:
        merged_data[gene] = stat[gene]
        merged_data[gene]["gc3"] = fold4[gene]["gc3"]
        merged_data[gene]["polNS"] = pnps[gene]["pns"]
        merged_data[gene]["polS"] = pnps[gene]["ps"]
        
    return merged_data
    

def conv_to_np(data):
    header = list(data[next(iter(data))].keys())
    data_array = np.zeros((len(data), int(len(header))))
    gene_order = [0]* len(data)

    for x, gene in enumerate(data):
        gene_order[x] = gene
        data_array[x,:len(header)] = [data[gene][y] for y in header]
    return {"gene": gene_order, "header": header, "metrics": data_array.astype(float)}


    
 



sumstat = parsing_input(sumstat_file)
fourfold = parsing_input(fourfold_file)
pnps = parsing_input(PnPs_file)

merged_dataset = merging_dataset(sumstat, fourfold, pnps)
final = conv_to_np(merged_dataset)

###outputing files
out = open(output_file, "w")
out.write("name" + "\t" + "\t".join(final["header"]) + "\n")
final_str = list(zip(final["gene"], ["\t".join(list(map(str,x))) for x in final["metrics"]]))
out.write("\n".join(["\t".join(x) for x in final_str]))




