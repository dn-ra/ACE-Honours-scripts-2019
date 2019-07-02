#!/srv/sw//python3/3.5.0/bin/python3


# -*- coding: utf-8 -*-
"""
Spyder Editor

Author: Daniel Rawlinson, Australian Centre for Ecogenomics
Email: daniel.rawlinson@uqconnect.edu.au

This script is made to collate the hits from a blast result and determine the quantities of plasmid vs chromosomal origin.
Elements stolen from Kai Blin: https://github.com/kblin/ncbi-acc-download
"""


from argparse import ArgumentParser
import requests
import xml.etree.ElementTree as ET
import errors
import csv
import sys

try:
    from httplib import IncompleteRead
except ImportError:
    from http.client import IncompleteRead

#stolen from Kai Blin: https://github.com/kblin/ncbi-acc-download
def get_stream(url, params):
    """Get the actual streamed request from NCBI."""
    try:
        r = requests.get(url, params=params, stream=True)
    except (requests.exceptions.RequestException, IncompleteRead) as e:
        print("Failed to download {!r} from NCBI".format(params['id']), file=sys.stderr)
        raise errors.DownloadError(str(e))

    if r.status_code != requests.codes.ok:
        print("Failed to download file with id {} from NCBI".format(params['id']), file=sys.stderr)
        raise errors.DownloadError("Download failed with return code: {}".format(r.status_code))
        
    return r


##os.join?





"""Command line handling."""
parser = ArgumentParser()
##input must be in form of csv IDs
##this currently exists for list of hits from each putative plasmid contig
parser.add_argument('ids', nargs=1, metavar='Blast output in csv format')
## for selection of database depending on the input sequence alphabet
parser.add_argument('-m', '--molecule', choices=["nucleotide", "protein"] ,help="Sequence type") #necessary? will proteins hit to whole genomes?
#to create and select directory for xml outputs
parser.add_argument('-d', '--dir', default='.', help='directory name for storage of outfile')
parser.add_argument('-i', '--infmt', choices=["blast", "diamond"], help='blast or diamond input filetype')

opts=parser.parse_args()

print('options invoked: \n\t input file: {} \n\t molecule tyep: {} \n\t output directory: {} \n\t input file type {}:'.format(opts.ids[0], opts.molecule, opts.dir, opts.infmt))

#make sure file writes firsts
if '/' in opts.ids[0]:
    filename = opts.ids[0].split('/')[-1]
else:
    filename = opts.ids[0]

if opts.dir:
    filename = "{}/{}_hitsummary.csv".format(opts.dir, filename)
else:
    filename = "{}_hitsummary.csv".format(filename)
open(filename, 'w+').close()

ENTREZ_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'
params = {'db': opts.molecule, 'version':2.0}

#extract blast hit IDs from plasmid prediction results

hits ={}
print('reading in blast output')

if opts.infmt == 'diamond':
    delim = '\t'
else:
    delim = ','

with open(opts.ids[0], 'r', encoding='utf-8') as f:
    readfile = csv.reader(f, delimiter = delim)
    for line in readfile:
        if line: #for some reason the reader object has empty lines read in?
            hits.setdefault(line[0], [])
            hits[line[0]].append(line[1])
        

#begin iteration over each contigs blast hits
genometypedict = {}            
            
contig_dict = {}
negative_taxids = [] #to get details of hits that aren't tagged as plasmid. Maybe they are viruses?


for key in hits.keys():
    ids = [[]] #id list for this contig, divided into separate lists in case of long url
    i=0 #index for storage of ids
    plascnt = 0 #start count for number of plasmid hits from this contig
    totalhits = 0 #count valid hits from each contig
    l=0 #length of urls
    for hit in hits[key]:
        l = l + len(hit) + 1 #+1 for length comma when joined
        if  l >= 2000- len(ENTREZ_URL): #approxomate url limit of 2000
            i+=1
            l=0
            ids.append([])
        ids[i].append(hit)
           
    
#retrive data
#for each contig BLAST results (which constitutes any multitude of hits)
    print('getting blast hits for contig {}'.format(key))

    for idlist in ids:
        params['id'] = ','.join(idlist)
    
    #Actual DocSum retrieval step
        with get_stream(ENTREZ_URL, params) as r:
       
    
            #TODO - set up to find genome tag with indexing rather than iterating
            #write records to element tree            
            contig_hits = ET.ElementTree(ET.fromstring(r.text))
            root = contig_hits.getroot()
            elems = [[elem for elem in child] for child in root] #take individual entries
            for child in elems[0]:
                if child.find('Genome') != None: #if entry not suppressed
                    temptag = child.find('Genome').text
                    if temptag != None:
                        totalhits+=1
                        if temptag in genometypedict:
                            genometypedict[temptag] +=1
                        else:
                            genometypedict[temptag] = 1
                    
                    if temptag == 'plasmid':
                        plascnt +=1 #count hits to plasmid
                ##else if not plasmid, get taxid
                    else:
                        negative_taxids.append(child.find('TaxId').text)
    #output contig-hits ratio dictionary            
    if totalhits == 0: #avoid div-by-0 error
        totalhits = .01
    
    contig_dict[key] = totalhits, plascnt, plascnt/totalhits
    

#Get Taxonomy information for non-plasmid hits
print('getting taxonomy information for non-plasmid hits')
params['db'] = 'taxonomy'
taxtag = {}
ranktag = {}

taxids = [[]] #id list for this contig, divided into separate lists in case of long url
i=0 #index for storage of ids
l=0 #length of urls
for tax in negative_taxids:
    l = l + len(tax) + 1 #+1 for length comma when joined
    if  l >= 2000- len(ENTREZ_URL): #approxomate url limit of 2000
        i+=1
        l=0
        taxids.append([])
    taxids[i].append(tax)
        
for lists in taxids:
    params['id'] = ','.join(lists)
    with get_stream(ENTREZ_URL, params) as r:
            tree = ET.ElementTree(ET.fromstring(r.text))
            root = tree.getroot()
            elems = [[elem for elem in child] for child in root] #take individual entries
            for child in elems[0]:
                if child.find('GenBankDivision') != None:
                    temptag = child.find('GenBankDivision').text
                    if temptag in taxtag:
                        taxtag[temptag] += 1
                    else:
                        taxtag[temptag] = 1
                if child.find('Rank') != None:
                    temptag = child.find('Rank').text
                    if temptag != None:
                        if temptag in ranktag:
                            ranktag[temptag] +=1
                        else:
                            ranktag[temptag] = 1
                            
    
#output plasmid hits
#TODO - include p-vales and total hit counts in output?

print('writing output files')

#sorteddict = sorted(contig_dict.items(), key=operator.itemgetter(0))
with open(filename, 'w') as f:
    w = csv.writer(f, delimiter = '\t')
    w.writerow(['ContigID', 'TotalHits', 'Plasmid Hit Count', 'Plasmid/Hit ratio'])
    for key, value in contig_dict.items():
        w.writerow([key, value[0], value[1], value[2]])

#output non-plasmid taxonomy ids
if opts.dir:
    filename = "{}/{}_nonplasmidtax.txt".format(opts.dir, filename)
else:
    filename = "{}_nonplasmidtax.txt".format(filename)
with open(filename, 'w') as f:
    w = csv.writer(f, delimiter = '\t')
    w.writerow(['Division', 'Count'])
    for div, count in taxtag.items():
        w.writerow([div, count])
    w.writerow([])
    w.writerow(['Rank', 'Count'])
    for rank, count in ranktag.items():
        w.writerow([rank, count])

print(taxtag)
print(ranktag)

print(genometypedict)
    

