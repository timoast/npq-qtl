#! /usr/bin/env python

import os
from glob import glob


for key, value in data.items():
    os.makedirs("../ProcessedData/" + key)
    fastqs = glob("./" + data[key][1] + "*/*.fastq.gz")
    for i in fastqs:
        trailing = i.split("/")[-1]
        os.rename(i, str("../ProcessedData/" + key + "/" + trailing))

for i in glob("./SRR*"):
    if not os.listdir(i):
        os.rmdir(i)
