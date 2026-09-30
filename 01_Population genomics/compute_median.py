#!/bin/python3

"""
Compute median based on count data to avoid storing large arrays to memory
python3 "input file" "column number" "id column" "output file" "quantile" (default = 0.5)
"""

import sys
from collections import defaultdict


#Parsing user input
infilepath = sys.argv[1]
datacolumn = int(sys.argv[2])
idcolumn = int(sys.argv[3])
output = sys.argv[4]
q = float(sys.argv[5]) if len(sys.argv) >= 6 else 0.5

def parse_median(count_dict, total,quantile = 0.5):
    if len(count_dict) == 1 :
        return list(count_dict)[0]

    cumtot = 0
    sorted_keys = sorted(count_dict)
    
    for y, i in enumerate(sorted_keys):
        cumtot += count_dict[i] / total

        if abs(quantile - cumtot) <= 1/total:
            res = (i + sorted_keys[y + 1])/2
            return str(int(res)) if int(res) == res else str(res)


        if cumtot > quantile:
         return str(i)


counts = defaultdict(lambda: defaultdict(int))
total_pos = defaultdict(int)
with open(infilepath, "r") as infile:
    for line in infile:
        l = line.strip().split("\t")
        counts[l[idcolumn]][int(l[datacolumn])] += 1
        total_pos[l[idcolumn]] += 1

results = {x: parse_median(counts[x], total_pos[x], q) for x in counts}
with open(output, "w") as out:
    out.write("\n".join([f"{x}\t{results[x]}" for x in results]))
    
