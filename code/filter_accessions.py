#! /usr/bin/env python

keep_accessions = [line.rsplit()[0] for line in open("../RawData/accessions.txt", "r")]

data = {}
with open("../RawData/A_thaliana_master_accession_list.csv", 'r') as infile:
    next(infile)
    for line in infile:
        line = line.split(',')
        id = line[0]
        acc = line[1]
        if acc in keep_accessions:
            data[id] = acc

with open("../RawData/1001_genomes_sra_runinfo.csv", 'r') as infile, open('../RawData/accession_sra.tsv', 'w+') as outfile:
    next(infile)
    for line in infile:
        line = line.split(',')
        sra = line[0]
        id = line[29]
        if id in data.keys():
            outfile.write("\t".join([sra, id, data[id]]) + "\n")
