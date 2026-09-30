#! /usr/bin/env python

from __future__ import print_function

from collections import defaultdict
import gzip
import optparse
from time import time
import sys
START = time()

# Function to find information in position in FORMAT field
def findposinformatstr(info_str, str2_id):
    returnpos = 0
    info = info_str.split(":")
    for i in info:
        if i == str2_id:
            break
        returnpos = returnpos + 1
    return returnpos


def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

parser = optparse.OptionParser()
parser.add_option('-f', '--vcf_file', action="store", dest="namefp", type="string")


# Unpack user parameters
(options, args) = parser.parse_args()

# Open input vcf
vcfin = myopen(options.namefp)

# MAIN
# Init samples / chromosome names
nameSeq = []
nameInd = []
chromSize = dict()  # **NEW: dictionary that will contain contig length

# init position loop / memory of preceding lines
oldchrom = ""
countSNP = 0
NumberOfSeq = ""
sites = []
pos = 0
oldpos = 0

for line in vcfin:
    line = line.rstrip()
    # ***NEW: Parsing chromosome length to fill up potential gap at end of chromosome [IE:BCFTOOLS]
    if line.startswith("##contig"):
        chromSize[line[13:line.find(",")]] = int(line.split("=")[-1].strip(">"))

    # count number of individuals and store names in array
    if line.startswith("#CHROM") and len(nameInd) == 0:
        arrline = line.split()

        for i in arrline[9:]:
            indseq = i
            nameSeq.append(indseq)
            nameInd.append(i)

        NumberOfSeq = len(nameSeq)
        print("There is " + str(NumberOfSeq) + " individuals")
        continue

    if line.startswith("#"):
        continue

    # init infos variables
    arrline = line.split()
    chromosome = arrline[0]
    pos = int(arrline[1])
    REF = str(arrline[3])
    ALT = str(arrline[4])

    # init chromosome name 
    if oldchrom == "":
        oldchrom = chromosome


    # ***NEW: for indels, bcftools sometimes report several time the same position
    if pos == oldpos and chromosome == oldchrom:
        print(f"Multiple entry for one position, first entry kept: {chromosome}-{pos}")
        continue
    # ***End of New
    
    if chromosome != oldchrom:  # HANDLES NEW CHROMOSOMES
        # **NEW: when switch to a new chromosome, add "N" at the end
        # if last parsed position is not equal to chromosome length

        if oldpos < chromSize[oldchrom]:
            numberSite = int(chromSize[oldchrom]) - oldpos
            [sites.extend(["N"] * len(nameSeq)) for _ in range(numberSite)]
       # **ENDofNEW
        

        fasta = open(oldchrom + ".fst", "w")
        for i, ind in enumerate(nameSeq):
            fasta.write(f">{ind}\n{''.join(sites[i::NumberOfSeq])}\n")
        
        oldchrom = chromosome
        oldpos = 0  # ***NEW: reset oldpos to 0 to handle new chromosome not begining at 0
        sites = []
        seq = defaultdict()
    
    # MODIF PROPOSER PAR EMERIC START UPGRADED BY FLORENT###############################
    # This was at the end before
    if (oldpos + 1) < pos:  # if current site is not exactly one position after oldpos site put "N"
        numberSite = int(pos - oldpos)
        [sites.extend(["N"] * len(nameSeq)) for _ in range(1, numberSite)]
    # MODIF PROPOSER PAR EMERIC STOP###############################
    

    # print(time() - START)
    covPos = findposinformatstr(arrline[8], "DP")  # find coverage position in INFO field
    gqPos = findposinformatstr(arrline[8], "GQ")

    if ALT == ".":
        for i in range(len(nameInd)):
            sites.extend([REF])
        oldpos = pos
        continue

    for i in range(0, len(nameInd)):
        arrInd = arrline[i + 9].split(":")
        sites.extend(arrInd[0].replace("0", REF).replace("1", ALT).replace(".", "N"))
    oldpos = pos
    # print(time() - START)


# ***NEW: when switch to a new chromosome, add "N" at the end if last parsed position is not equal to chromosome length
if oldpos < chromSize[oldchrom]:
    numberSite = int(chromSize[oldchrom] - oldpos)
    [sites.extend(["N"] * len(nameInd)) for _ in range(0, numberSite)]
    print("B")
# ENDofNEW


print("Print  " + oldchrom)
fasta = open(oldchrom + ".fst", "w")
for i, ind in enumerate(nameSeq):
    fasta.write(f">{ind}\n{''.join(sites[i::NumberOfSeq])}\n")

vcfin.close()
END = time()

print(END-START)
