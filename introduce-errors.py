#!/usr/bin/python

"""
Created on Fri May 24 13:22:51 2019

A script to introduce errors into a sequence without fragmenting and reassembling
Primary purpose is to test nucmer for aligning a reference sequence back to itself after introducint errors

@author: Daniel Rawlinson
@email: daniel.rawlinson@uqconnect.edu.au
"""

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from random import randint

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('sequence', help="Fasta file with single sequence to be edited")
parser.add_argument('--errorpercent', '-e', type=float, help="percentage of errors to introduce to sequence")

opts =parser.parse_args()

##parse arguments

##read in sequence

##extract sequence string

##apply error model

##reenter string into 

record = SeqIO.read(opts.sequence, 'fasta')
seqlist = list(record.seq)
letters = ['G','C','T','A']

if opts.errorpercent > 1: opts.errorpercent /= 100

bases = len(seqlist)
errornum = int(opts.errorpercent * bases)

print('seq length is {}'.format(bases))

indx = []
print('editing {} bases'.format(errornum))
for num in range(errornum):
    indx.append(randint(0, bases))

for elem in indx:
    i = letters.index(seqlist[elem]) #TODO - i think there is an off-by-1 error in here
    i+= randint(1,3)
    if i>3:
        i -= 4
    seqlist[elem] = letters[i]
    

recorderror = SeqRecord(Seq(''.join(seqlist)))

recorderror.id = record.id + '_werrors'
recorderror.name = recorderror.id
recorderror.description = record.id + ' with error ratio of {} added'.format(opts.errorpercent)

SeqIO.write(recorderror, recorderror.id, 'fasta')

print('Done!')
    
    

    


