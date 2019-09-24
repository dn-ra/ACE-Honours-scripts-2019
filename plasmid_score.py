#!/srv/sw//python3/3.5.0/bin/python3

"""
Created on Mon Sep  9 12:50:18 2019

@author: dan_r
"""

import csv
import argparse

parser = argparse.ArgumentParser(description='Score contigs on plasmid likelihood')
parser.add_argument('-p', help='prokka .tbl file')
parser.add_argument('-out', help='prefix to save files to')
args = parser.parse_args()

#constants
score_file = '/srv/projects/abisko/dan/repeatm_tests/all_assemblies/cluster_test_unionfind/pfamvscore'
#read in score file
pfam_score = {}
f = open(score_file, 'r')
r = csv.reader(f)
for line in r:
    pfam_score[line[0]] = float(line[1])



#read in prokka file
cdscluster = {}
with open(args.p, 'r') as f:
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

filename = "{}.clumped".format(args.out)

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



#process dictionary with plasmid-pfams data
filename="{}.scored".format(args.out)

contig_scores = {}
no_hits_count = 0


#these are for the enrichment data
plasmid_gene_count = 0
total_gene_count = 0
contig_enrichment = []


for contig, orfs in cdscluster.items():
    score = 1
    trna_count = 0
    rrna_count = 0
    if orfs:
        plasmid_genes = 0
        total_genes = 0
        for info in orfs.values():
            #remove suffix from pfam accession
            pfam = info[1].split('.')[0]
            if pfam in pfam_score:
                score*=pfam_score[pfam]
                if pfam_score[pfam] > 10:
                    plasmid_genes+=1
                    plasmid_gene_count +=1
                else:
                    total_gene_count+=1
                    total_genes+=1
            elif info[0] == 'tRNA':
                trna_count+=1
            elif info[0] =='rRNA':
                rrna_count+=1
        
        contig_scores[contig] = [score, trna_count, rrna_count]
        contig_enrichment.append(plasmid_genes/total_genes) 
    else:
        no_hits_count+=1

f = open(filename, 'w')
w = csv.writer(f, delimiter = '\t')
w.writerow(["CONTIG ID", "PLASMID SCORE", "tRNAs", "rRNAs"])
for contig, scores in contig_scores.items():
    s1, s2,s3 = scores
    w.writerow([contig, s1,s2,s3])
    
    
#write plasmid enrichment value
filename = '{}.enriched'.format(args.out)

f = open(filename, 'w')
f.write('#{} plasmid genes counted above score of ten\n'.format(plasmid_gene_count))
f.write('#{} total non-hypothetical genes counted\n'.format(total_gene_count))
f.write('#Enrichment score: {}\n'.format(plasmid_gene_count/total_gene_count))
f.write(",".join(contig_enrichment))
        
f.close()

print('plasmid search performed. Scores stored in {}. {} contigs found with no identifiable orfs.'.format(filename, no_hits_count))
    
    
    

