'''This from RepeatM Clusterer module - author wwood.
To be used here for testing of single_contig clustering from delta match files'''

'''
    imports
        '''
        

import re
from logging import warning #also import exceptions?
from operator import itemgetter


#pattern to determine if contig name is in spades format
spadespattern = re.compile(r'NODE_[0-9]*_', re.UNICODE)

'''--------------------Begin class definition------------------------------'''
class Contig_Cluster(object):
    def __init__(self, node_list, matches):
        if isinstance(node_list, str):
            raise RuntimeError('input to Contig_Cluster class must be list. String has been entered.')
        self.nodes = node_list
        self.size = len(self.nodes)
        self.av_cov = None
        self.av_length = None
        self.matches = set(matches) 
         #has_spades necessary?       
        if self.has_spades() == True:
            length_total = 0
            cov_total = 0
            spades_seqs = 0
            for node in self.nodes:
                if spadespattern.match(node):
                    nodeinfo = node.split('_')
                    length_total += int(nodeinfo[-3])
                    cov_total += float(nodeinfo[-1])
                    spades_seqs +=1
                else:
                    warning('Some contigs in this cluster were not generated by Spades. Non-Spades contig names not supported yet')
            self.av_cov = cov_total / spades_seqs
            self.av_length = length_total / spades_seqs
        else:
            warning('Non-Spades contig names not supported yet')
    #necessary?
    def has_spades(self):
        #TODO - stop interpreting single contigs as list
        spades = False
        for c in self.nodes:
            if spadespattern.match(c):
                spades = True
                break
        return spades
    #TODO -
    def plas_label(self):
        return None
    
    
    def retrieve_seqs(self):
        #all I need are nodes and location of source fasta files? use sequence library in repeatM
        #pop out assembly number from start of contig? use as input the source fastafiles?
        nodestring = ''
        assemblystring = ''
        for n in self.nodes:
            nodesplit = n.split("__")
            nodestring+=(">"+nodesplit.pop()+"\n") #second entry into nodelist
            assemblystring+=(nodesplit+"\n") #remainder (1st entry) into assembly list
        
        #awk (retrieve sequence) $contigID $assembly_file
        nodefind_awk = '''-v name=$node 'BEGIN{RS=">";FS="\n"}NR>1{if ($1~/name/) print ">"$0}'''
        CMD = 'while read -r node <&1 && read -r assembly <&2; do awk {} $assembly; done 1<(printf "{}") 2<(printf "{}")'.format(nodefind_awk, nodestring, assemblystring)

        ##TODO - run subprocess
        return None
    
    def inherit_matches(self, match_dict): #apply to cluster, but inherit from dict_threshold output
        #for n in self.nodes:
           #write function here or call another function? 
         return None  
    #TODO -
    def label_cluster(self):
        '''label cluster as linear or circular based on alignment evidence
        Evidence includes: 
        - Does the sequence consistenly start and end with the same ORF? (don't need to know the proteins: Use OrfM)
        - Do regions align globally or differentially?
        - do endings overlap?
        '''
    #OrfM is already a subprocess in clusterer module
        return None



'''---------------------End class definition-------------------------------'''


def sort_clusters(cluster_list): #or maybe a dictionary instead?
    '''sort clusters by: 
    1. N in cluster, 2. length of N, 3. coverage of N
    Don't get it just from the names. That's SPAdes format but some people won't use spades
    Get it from the data directly. But megahit doesn't include any of this informaiton in the contig name, 
    and coverage would need reads mapped'''
    #TODO - sort out lambda functions
    sortbysize = lambda c: (c.size is not None, c.size)
    sortbylength = lambda c: (c.av_length is not None, c.av_length)
    sortbycov = lambda c: (c.av_cov is not None, c.av_cov)
    
    
    #main sort function
    return cluster_list.sort(reverse=True, key= lambda c: (sortbysize(c), sortbylength(c), sortbycov(c)))

def cluster_nucmer_matches(sig_matches): #sig_matches comes out of delta_parse.collate_sig_matches()
    sig_match_dict = {}
    match_links = []
    for m in sig_matches:
        link = m.seqs
        try:
            sig_match_dict[link[0]].append(m) 
            sig_match_dict[link[1]].append(m)
            
        except KeyError:
            sig_match_dict[link[0]] = [m] 
            sig_match_dict[link[1]] = [m]
            
        match_links.append(link)
    
    cluster_list = single_linkage_cluster(match_links) #main clustering step
    
    for cluster in cluster_list:
        nucmer_match_in_cluster = []
        for node in cluster:
            nucmer_match_in_cluster += set(sig_match_dict[node])
        Contig_Cluster(cluster, nucmer_match_in_cluster)     
        
    
    
    ''' Index here or later?
    sig_match_dict = {}
    for match in sig_matches:
        sig_match_dict[match.seqs[0]] = match
      '''  

'''!!! Don't edit this. This is what will be in the clusterer module of repeatm!!!'''
def single_linkage_cluster(links): #dictionary or list input? #originally a class function. changed here to just be 'links' as input for testing. Refer to original RepeatM module for original inputs.
    '''Not sure this is the fastest, but eh. Return a list of lists, where each
    list is a cluster.'''
    clusters = {}
    next_group_number = 1 # counter of cluster numbers. goes up by one as each new cluster grouping is created
    for link in links: #for each link generated by mash and such
        one = link[0]
        two = link[1] #first and second elements of the match
        if one in clusters and two in clusters: #see if both of these are already clustered somewehre. bring those clusters together
            old_one_number = clusters[one] # ._._. set 'old_one' as cluster[one], set cluster[one] as cluster[two]
            clusters[one] = clusters[two] #now both one and two have same cluster, old cluster of one is stored in tmp variable
            names_to_fix = [] #??
            for name, number in clusters.items():
                if number == old_one_number: 
                    names_to_fix.append(name) #get all names that are already in that old cluster and append them to a list
            for n in names_to_fix:
                clusters[n] = clusters[two] #set cluster number for all the names_to_fix sequences to be the same as clusters[two] (which is now also the same as clusters[one]
        elif one in clusters:
            clusters[two] = clusters[one] # of one already clustered, cluster 2 with it
        elif two in clusters:
            clusters[one] = clusters[two] # if two already clustered, cluster 1 with it
        else:
            clusters[one] = next_group_number
            clusters[two] = next_group_number
            next_group_number += 1 # if not clustered, set both in new cluster, increase cluster counter by one for next new cluster
    reverse = {} #dictionary of reversed name: cluster dict
    for name, number in clusters.items():
        if number not in reverse:
            reverse[number] = [] #add that cluster number into dict if not there yet
        reverse[number].append(name) #add id into that cluster dict
    return list(reverse.values()) #output as list of lists I GET IT!
    