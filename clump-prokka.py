# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 12:50:18 2019

@author: dan_r
"""

import csv


file = 'PROKKA_09062019.tbl'
fileprefix = 'Prokka'

#read in prokka file
cdscluster = {}
with open(file, 'r') as f:
    readfile = csv.reader(f, delimiter='\n')
    for line in readfile:
        if line[0].startswith('>'):
            currentkey = line[0].replace(">Feature ", "")
            cdscluster.setdefault(currentkey, [])
        if line[0].startswith('\t\t\tlocus'):
            cdscluster[currentkey].append(line[0].split('\t')[-1])

#write clustering to file

filename = "{}.clumped".format(fileprefix)

with open(filename, 'w') as f:
    w = csv.writer(f, delimiter = '\t')
    #header here
    w.writerow(["CONTIG ID", "ORF ID", "FIRST TAXID", "GENOME LABEL", "LABELLED ID", "ANNOTATION"])
    for key, value in cdscluster.items():
        w.writerow([key])
        if value:
            for elem in value:
                w.writerow(['\t' + elem])

        else:
            w.writerow(['\t' + '--No hits--'])  # (3)
