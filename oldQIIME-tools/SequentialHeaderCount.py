#!/usr/bin/python3
# A program for modifying fasta headers
# USAGE: ./changingHeaders_devon.py
# Author: Taruna Aggarwal
# Affiliation: University of New Hampshire, Durham, NH, USA
# Date: 1/17/2016
# Purpose is to replace parts of headers and add consecutive numbers to each sample ID

import sys
from collections import Counter

inFile = open(sys.argv[1], "r")
outFile = open(sys.argv[2], "w")
counter = Counter()

for line in inFile:
    line = line.rstrip()
    if line[0] == ">":
        sample_id = line.split('\t')[0]
        counter[sample_id] += 1
        outFile.write("{0}_{1}\n".format(sample_id, str(counter[sample_id])))
    else:
        outFile.write("{0}\n".format(line))
