#!/srv/sw//python3/3.5.0/bin/python3


# -*- coding: utf-8 -*-


"""
TODOS 
    - Save hit length in final output
    - distinguish between ORF that have no hits and ORF hits that don't have tags
    - rename final outputs
    - Summarise ORFs for each contig
    - change so that  prokka isn't needed
    """


"""
Spyder Editor

Author: Daniel Rawlinson, Australian Centre for Ecogenomics
Email: daniel.rawlinson@uqconnect.edu.au

This script is made to collate the hits from a blast result and determine the quantities of plasmid vs chromosomal origin.
Elements stolen from Kai Blin: https://github.com/kblin/ncbi-acc-download
"""

"""library imports"""

from argparse import ArgumentParser
from collections import OrderedDict
from time import sleep
import requests
import xml.etree.ElementTree as ET
import errors
import csv
import sys
import os

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


"""Command line handling."""
parser = ArgumentParser()
##input must be in form of blast or diamond outupt in default tsv format
##this currently exists for list of hits from ORFs on predicted plasmid contigs
parser.add_argument('ids', nargs=1, metavar='Blast or diamond output in default tsv format')
## for selection of database depending on the input sequence alphabet
parser.add_argument('-m', '--molecule', choices=["nucleotide", "protein"] ,help="Sequence type") #necessary? will proteins hit to whole genomes?
#to create and select directory for xml outputs
parser.add_argument('-d', '--dir', default='.', help='directory name for storage of outfiles')
parser.add_argument('-c', '--cluster', help= 'must be .tbl file from prokka output. Used for clustering CDS around parent contig')

opts=parser.parse_args()

"""feedback to user on options"""
print('options invoked: \n\t input file: {} \n\t molecule tyep: {} \n\t output directory: {} \n\t clustering? {}'.format(opts.ids[0], opts.molecule, opts.dir, 'Yes' if opts.cluster else 'No'))
if opts.cluster:
    print('\t\t clustering with {}'.format(opts.cluster))

#TODO - tidy this all up
#make sure file writes firsts
fileprefix = opts.ids[0].split('/')[-1]
fileprefix = os.path.splitext(fileprefix)[0]

if opts.dir:
    if not os.path.isdir(opts.dir):
        os.mkdir(opts.dir)
    filename = "{}/{}_cdstags.csv".format(opts.dir, fileprefix)
else:
    filename = "{}_cdstags.csv".format(fileprefix)

with open(filename, 'w') as f:
    w = csv.writer(f, delimiter = '\t')
    w.writerow(['CDS ID', 'Genome Tag', 'Accession #', 'Annotation'])

ENTREZ_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'
params = {'db': opts.molecule, 'version':2.0}
tax_params = {'db': 'taxonomy', 'version': 2.0}
#extract blast hit IDs from plasmid prediction results

hits =OrderedDict()
print('reading in blast outupt')


delim = '\t' #tab separated input file
with open(opts.ids[0], 'r', encoding='utf-8') as f:
    readfile = csv.reader(f, delimiter = delim)
    for line in readfile:
        if line: #for some reason the reader object has empty lines read in?
            hits.setdefault(line[0], [])
            hits[line[0]].append(line[1])
        

#begin iteration over each ORFs blast hits
orf_dict = {}
tax_dict = {}

for key in hits.keys(): #for each ORF
    
    write = False #counter to keep in loop until a valid hit tag is found for the ORF
    tax = False #switcher to get taxid of first hit
    ids = [[]] #id list for this ORF, divided into separate lists in case of long url
    i=0 #index for storage of ids
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
    print('retrieving hit info for ORF {}'.format(key))
    i=1
    for idlist in ids:
        print('running idlist {} of {}'.format(i, len(ids)))
        params['id'] = ','.join(idlist)
    
    #Actual DocSum retrieval step
        
        for attempt in range(5):
            repeat = None #Set to break loop at end if successful

            try:
                with get_stream(ENTREZ_URL, params) as r:
                    #write records to element tree
                    orf_hits = ET.ElementTree(ET.fromstring(r.text))
                    root = orf_hits.getroot()
                    elems = [[elem for elem in child] for child in root] #take individual entries
                    for child in elems[0]:
                        if tax == False:
                            temptax = child.find('TaxId').text
                            if temptax:
                                tax_dict[key] = temptax
                                tax = True
                                print(key, tax)
                            else:
                                print(key, 'no taxtag')
                        if child.find('Genome') != None: #if entry not suppressed
                            temptag = child.find('Genome').text
                            if temptag != None:
                                print('writing', temptag)
                                hit = child.find('Caption').text
                                annot = child.find('Title').text
                                orf_dict[key] = temptag, hit, annot
                                write = True
                                break #return to idlist loop
                            else:
                                orf_dict[key] = '--No labels--'
                                #return to idlist loop
                                i+=1
                        else:
                            pass
            #continue if temptag is still none and there are more idlists to go through

                    
            #TODO include in record: hit name, hit index (ie. the hit number for that CDS)
            
               
        
            except ConnectionError:
                print('Connection Error. Sleeping...')
                sleep(2)
                repeat = True
                pass
                #retry          
            except errors.DownloadError:
                print('Download Error. Retrying...')
                repeat = True
                pass
                #retry
            
            
            if repeat == True:
                continue #go through repeat loop again
            else:
                
        
                    
                break #no errors - do not repeat
        
    #write record to file. Is this necessary?
    
    with open(filename, 'a') as f:
        w = csv.writer(f, delimiter = '\t')
        if write == True:
            w.writerow([key, temptag, hit, annot])
        else:
            w.writerow([key])

print('hit tags made and written to file')

print(tax_dict)

"""cluster CDS tags together around contig of origin"""

if opts.cluster:
    print('clustering...')
    #read in PROKKA tbl file
    cdscluster = {}
    with open(opts.cluster, 'r') as f:
        readfile = csv.reader(f, delimiter='\n')
        for line in readfile:
            if line[0].startswith('>'):
                currentkey = line[0].replace(">Feature ", "")
                cdscluster.setdefault(currentkey, [])
            if line[0].startswith('\t\t\tlocus'):
                cdscluster[currentkey].append(line[0].split('\t')[-1])
    
    #write clustering to file
    #3 outcomes to include: 1. No ORFs on the contig, 2. No labels in the recovered hits, 3. No hits from DMND at all
    with open(fileprefix +'.clustered.tsv', 'a') as f:
        w = csv.writer(f, delimiter = '\t')
        #header here
        w.writerow(["CONTIG ID", "ORF ID", "FIRST TAXID", "GENOME LABEL", "LABELLED ID", "ANNOTATION"])
        for key, value in cdscluster.items():
            division_dict = {}
            w.writerow([key])
            if not value:
                w.writerow(['\tNo ORFs on this contig']) #if there are no ORFs for that contig (1)
            else:
                for elem in value:
                    taxids = []
                    if elem in orf_dict:
                        cdstags = orf_dict[elem]
                        taxtag = tax_dict[elem]
                        taxids.append(taxtag)
                        if isinstance(cdstags, str):
                            w.writerow(['\t'+ elem, taxtag, cdstags]) #This should say --No labels-- (2)
                        else:
                            w.writerow(['\t' + elem, taxtag, cdstags[0], cdstags[1], cdstags[2]])

                    else:
                        w.writerow(['\t' + elem, '--No hits--'])  # (3)
            
            #retrieve taxid information and write
            tax_params['id'] = ",".join(taxids)
#            division_dict = {}
            with get_stream(ENTREZ_URL, tax_params) as t:
                tax_info = ET.ElementTree(ET.fromstring(t.text))
                root = tax_info.getroot()
                elems = [[elem for elem in child] for child in root]
                for child in elems[0]:
                    if child.find('GenBankDivision') != None:
                        entry = child.find('GenBankDivision').text
                        try:
                            division_dict[entry] +=1
                        except KeyError:
                            division_dict[entry] = 1
            w.writerow(['Top hits taxonomy: ', division_dict])