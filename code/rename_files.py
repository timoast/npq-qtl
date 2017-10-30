#! /usr/bin/env python

import os
from glob import glob


data = {}
with open("../RawData/accession_sra.tsv", "r") as infile:
    for line in infile:
        line = line.rsplit()
        data[line[0]] = [line[1], line[2]]

for key, value in data.items():
    os.makedirs("../ProcessedData/" + key)
    fastqs = glob("./" + key + "*/*.fastq.gz")
    for i in fastqs:
        trailing = i.split("/")[-1]
        os.rename(i, str("../ProcessedData/" + key + "/" + trailing))

for i in glob("./SRR*"):
    if not os.listdir(i):
        os.rmdir(i)
