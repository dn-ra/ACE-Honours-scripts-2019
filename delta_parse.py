# -*- coding: utf-8 -*-
"""
A library for parsing and handling of information form nucmer's delta outfile. For use in RepeatM contig matching.

Created on Thu May 30 13:42:15 2019

@author: Daniel Rawlinson, Australian Centre for Ecogenomics (ACE)
@email: daniel.rawlinson@uqconnect.edu.au
"""


'''
imports?
    '''




'''------------------------begin class definition---------------------------'''
class Nucmer_Match(object):
    '''a format in which to store each separate nucmer sequence alignment'''
    seqs = None
    lengths = None
    hitstarts_1 = None
    hitstops_1 = None
    hitstarts_2 = None
    hitstops_2 = None
    numerrors = None
    name = None
    
    
    ''' self is the nucmer_match object, match is the seq names and seq lengths,
    match deets is a dictionary containing separate arrays of hit starts (x2), hit stops (x2),
    and numerrors'''
    def __init__(self, name, match, matchdeets):
        
        #TODO - need float here or integer? integer only works for python 3
        self.seqs =        [match[0],match[1]]
        self.lengths =     list(map(int, [match[2], match[3]]))
        self.hitstarts_1 = list(map(int, [entry[0] for entry in matchdeets]))
        self.hitstops_1 =  list(map(int, [entry[1] for entry in matchdeets]))
        self.hitstarts_2 = list(map(int, [entry[2] for entry in matchdeets]))
        self.hitstops_2 =  list(map(int, [entry[3] for entry in matchdeets]))
        self.numerrors =   list(map(int, [entry[4] for entry in matchdeets]))
        self.name = name
        
    def __len__(self):
        '''what happens when you call length of the class object'''
        return len(self.numerrors)
        
    def display(self):
        print('> '+(' '.join(self.seqs))+'\n')
        print('\t\tSeq1 region\tSeq2 region\tNumErrors')
        for i in range(len(self)):
            print('Match {}\t\t{}:{}\t\t{}:{}\t\t{}'.format(i+1, self.hitstarts_1[i], self.hitstops_1[i], self.hitstarts_2[i], self.hitstops_2[i], self.numerrors[i]))


    def gen_statistics(self):
        '''produce length ratio and bi-direcitonal alignment statistics to measure closeness of sequences''' 
        #TODO - rewrite so that length_ratio is always fraction
        #One variable in each converted to float so division works properly in Python 2
        length_ratio = float(self.lengths[0]) / self.lengths[1]
        
        #TODO - fix abs to put in correct place
        matchlength_1 = sum([float(self.hitstops_1[i]) - self.hitstarts_1[i] for i in range(len(self))])
        alignment_ratio_1 = abs(matchlength_1/self.lengths[0])
        
        matchlength_2 = sum([float(self.hitstops_2[i]) - self.hitstarts_2[i] for i in range(len(self))])
        alignment_ratio_2 = abs(matchlength_2/self.lengths[1])        
           
        
        return length_ratio, alignment_ratio_1, alignment_ratio_2
        pass


    
    
    #TODO - function for detecting and dealing with overlap situations

'''--------------------------end class definition---------------------------'''


#'''Temp file for testing'''
#file = 'C:\\Users\\dan_r\\Documents\\Honours_Data\\nucmertest\\testsplitvserrors\\splitvserrors_0.2_errors.delta'

''' Init arrays for alignment data of each match.
    Placed together, starts, stops, and errors will form a matrix with number of alignments equal to width
    '''

def deltaread(file):
    '''read through delta files'''
    with open(file, 'r') as f:
        recording = False #"recording" (boolean variable) is a switcher that will determine whether the line is being stored or not (ie. to ignore positions of indels)
        matchdeets = []
        name = None
        delta = {} #name of dictionary for storage of Nucmer_Match objects
        for line in f.readlines():
            if line.startswith('>'):
                if name: #skip flushing to dictionary if this is the first match record
                    #FIRST - flush previous hit to a nucmer_match object
                    delta[name] = Nucmer_Match(name, match, matchdeets)            
                    #remove alignment details of previous match
                    del matchdeets[:]
                #THEN - read in match details
                recording = True
                match = line.replace('>', '').split()
                #TODO - implement a shorter unique naming system
                name = '_'.join(match)
                
                continue
            elif line == '0\n':
                recording = True
                continue
            elif recording == True:
                matchdeets.append(line.split()[:-2])
                
                recording = False
    
    #flush last match to  nucmer_match object
    delta[name] = Nucmer_Match(name, match, matchdeets)

    return delta

def apply_threshold(deltadict, threshold = 0.97, outfile = None):
    '''input = dictionary output from deltaparse function
    apply blanket ratio threshold level and select matches to send to FastANI'''
    thresh_matches = []
    for key, value in deltadict.items():
        stats = value.gen_statistics()
        if all(stat >= threshold for stat in stats):
            thresh_matches.append(list(value.seqs))
    
    '''optionally write to file'''   
    if outfile:
        write_thresh_matches(thresh_matches, outfile)
        
    return thresh_matches


def write_thresh_matches(thresh_matches, outfile):
    '''write simple txt file containing sequence matches that pass threshold'''
    with open(outfile, 'w') as f:
        for elem in thresh_matches:
            f.write(elem[0] +' '+elem[1]+'\n')

#TODO - function to export in JSON format?