#! /usr/bin/env python

import os
import sys
from glob import glob


names = {}
with open(sys.argv[1], "r") as infile:
    for line in infile:
        line = line.rsplit()
        names[line[0]] = line[2]

files = glob("*_npq*")
for i in files:
    sra = i.split("_")[0]
    os.rename(i, names[sra] + i.strip(sra))
