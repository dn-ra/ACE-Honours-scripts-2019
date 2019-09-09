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
        if line[0].startswith('>'): #indicates beginning of new contig
            currentkey = line[0].replace(">Feature ", "")
            cdscluster.setdefault(currentkey, {})
        elif line[0].split()[0].isdigit(): #find CDS or tRNA tag
            if 'gene' not in line[0]:
                typ = line[0].split()[2]
        elif line[0].startswith('\t\t\tinference'):
            line = next(readfile) #go to next line to get inference details. if no second inference line, it's a hypothetical protein
        
        if line[0].startswith('\t\t\tinference'):
            pfam = line[0].split(':')[-1]
        elif line[0].startswith('\t\t\tlocus'):
            locus =line[0].split('\t')[-1]
        elif line[0].startswith('\t\t\tproduct'):
            product = line[0].split('\t')[-1]
            if product != 'hypothetical protein':
                cdscluster[currentkey][locus] = [typ, pfam, product]
#write clustering to file

filename = "{}.clumped".format(fileprefix)

with open(filename, 'w') as f:
    w = csv.writer(f, delimiter = '\t', quoting = csv.QUOTE_NONE, escapechar = '\\')
    #header here
    w.writerow(["CONTIG ID", "ORF ID", "TYPE", "PFAM", "ANNOTATION"])
    for key, value in cdscluster.items():
        w.writerow([key])
        if value:
            for locus, product in value.items():
                w.writerow(['\t' + locus + '\t' + '\t'.join([product[0], product[1], product[2]])], )
        else:
            w.writerow(['\t' + '--No hits--'])  
