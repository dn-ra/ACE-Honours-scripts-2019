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
    trna = False #init boolean
    for line in readfile:
        if line[0].startswith('>'):
            currentkey = line[0].replace(">Feature ", "")
            cdscluster.setdefault(currentkey, {})
        elif line[0].split()[0].isdigit():
            if 'tRNA' in line[0]:
                trna = True
                print('foundtrna')
        elif line[0].startswith('\t\t\tlocus'):
            locus =line[0].split('\t')[-1]
        elif line[0].startswith('\t\t\tproduct'):
            product = line[0].split('\t')[-1]
            cdscluster[currentkey][locus] = [product]
            cdscluster[currentkey][locus].append(str(trna))
            trna = False #reset boolean
#write clustering to file

filename = "{}.clumped".format(fileprefix)

with open(filename, 'w') as f:
    w = csv.writer(f, delimiter = '\t')
    #header here
    w.writerow(["CONTIG ID", "ORF ID", "ANNOTATION", 'is tRNA'])
    for key, value in cdscluster.items():
        w.writerow([key])
        if value:
            for locus, product in value.items():
                w.writerow(['\t' + locus + '\t' + '\t'.join([product[0], product[1]])], )
        else:
            w.writerow(['\t' + '--No hits--'])  
