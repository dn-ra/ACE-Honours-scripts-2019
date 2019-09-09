"""
A module for use in RepeatM. See: github.com/wwood/RepeatM

@author: Daniel Rawlinson, Australian Centre for Ecogenomics (ACE)
@email: daniel.rawlinson@uqconnect.edu.au
"""


'''
    imports
        '''
        
import os
import re
import subprocess
from logging import warning #also import exceptions?
import csv
import union_find_cluster



'''constants'''
#pattern to determine if contig name is in spades format
spadespattern = re.compile(r'.*NODE_[0-9]*_', re.UNICODE)

#set delimiter here?
#delim = '__'

#set sam_cmd here?

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
            self.av_cov = int(cov_total / spades_seqs)
            self.av_length = int(length_total / spades_seqs)
        else:
            warning('Non-Spades contig names not supported yet')
            
            
    def get_rep_seq(self):
        top_len = 0
        match_nodes = [match.seqs for match in self.matches]
        for n in self.nodes:
            match_num = len([m for m in match_nodes if n in m])
            if match_num > top_len:
                top_len = match_num
                top_n = n
        return top_n
    
    #necessary?
    def has_spades(self):
        #check if there is at least one spades contig in node_list
        spades = False
        for c in self.nodes:
            if spadespattern.match(c):
                spades = True
                break
        return spades
    
    def minimal_match(self):
        '''determine lowest passing statistics from matches in cluster
        In order of: length ratio, align a > b, align b > 1, ANIm'''
        low = [1,1,1,1]
        for m in self.matches:
            stats = m.gen_statistics()
            for i in range(len(stats)):
                if stats[i] < low[i]:
                    low[i] = stats[i]
            
        return low
            
            
            
    #TODO - for running plasflow or cbar?
    def plas_label(self):
        return None
    
    def split_names(self):
        '''splits the names of the nodes to redraw node-assembly links in dictionary format'''
        node_assembly_dict = {}
        for n in self.nodes:
            nodesplit = n.split("__")
            node_assembly_dict[">"+nodesplit[1]] = nodesplit[0] #remainder (1st entry) into assembly list
    
        return node_assembly_dict
    

    def retrieve_seqs(self, assembly_dir, repseq = False):
        #all I need are nodes and location of source fasta files?
        #pop out assembly number from start of contig? use as input the source fastafiles?
        #Beware leaked processes
        #can get assembly_dir from delta output?
        
        outdir = 'CLUSTER_size_{}_avlen_{}_avcov_{}'.format(self.size, self.av_length, self.av_cov)
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        
        if repseq == True:
            outfile = 'repseq_out.fa'
            assembly, node = self.get_rep_seq().split("__")
            nodes = {">"+node: assembly}
            
        else:
            outfile = 'cluster_out.fa'
            nodes = self.split_names()
                   
        seq_out = open(outdir + "/" + outfile, 'w')
        
        #parallelise eventually. That's why i've written it to dictionaries first
        for node, assembly in nodes.items():
            f = open("/".join([assembly_dir, assembly]))
            for i in f:
                if i.find(node)!=-1:
                    seq_out.write('>'+assembly+"__"+node.replace('>','')+'\n')
                    wholeseq =False #flag to tell me if I've taken the whole sequence yet
                    while wholeseq ==False:
                        line = next(f)
                        if line.startswith('>'):
                            wholeseq = True
                        else:
                            seq_out.write(line)
        
        
        seq_out.close()
        
        return None
    
#     previously proposed method
#       awk (retrieve sequence) $contigID $assembly_file
#        nodefind_awk = '''-v name=$node 'BEGIN{RS=">";FS="\n"}NR>1{if ($1~/name/) print ">"$0}\''''
#        CMD = 'while read -r node <&3 && read -r assembly <&4; do awk {} $assembly; done 3< <(printf "{}") 4< <(printf "{}")'.format(nodefind_awk, nodestring, assemblystring)
#        subprocess.call(['bash', '-c', CMD])
     
    
    
#what if I look for orientation of reads mapped from other assemblies? Shouldn't that be consistent too?
#have to call from within retreieve_seqs for tempfiles to be accessible
    def gen_minibam(node_assembly_dict, bam_location, outdir = 'minibam_out'): #pass self into this?
        ''' make this automated so that python can interpret the output itself?
        Will need to use Melody's script, plus something else I make to check for linearity'''
        sam_cmd = 'samtools view {1} {2} > %s/{1}_mini.bam' % (outdir) #edit this for correct command
        
        node_bam_dict = {}
        all_bams = [b for b in os.listdir(bam_location) if b.endswith('bam')]
        node_string = ''
        bam_string = ''
        for node, assembly in node_assembly_dict.items():
            bam_file = [b for b in all_bams if assembly.replace('fasta','') in b][0]
            node_bam_dict[node] = bam_file #not really necessary as I'm passing node-bam links straight into a string
            bam_string += bam_file+"\n"
            node_string += node+"\n"
            #removed input of temp_fna files
        #join into string separated by newline stuff to feed through to samtools is now in one long string with newline separators
        
    
        #TODO - load samtools and parallel into environment
        subprocess.call(['bash', '-c', 'parallel -a <(printf {}) -a <(printf {}) --link {}'.format(bam_string, node_string, sam_cmd)])
        
        return None
    #TODO -
    def label_cluster(self):
        ''' label cluster as linear or circular based on alignment evidence
        Evidence includes: 
        - Do regions align globally or differentially? #no of matches in Nucmer_Match
        - do endings overlap?
        -are there any regions that don't match?
        
        -classes of alignment. Applies to each match object
            + linear
                +perfect = ends of contigs align perfectly
                +imperfect = alignment does not extend perfectly to ends
                +complex = repeats or deletions
                +simple = no repeats or deletions
                    
            +circular
                +perfect
                +imperfect
                +complex
                +simple

Possibilites:
                            linear                |               circular
           |------------+-------------+-----------+-----------+-----------+-------
           | perfect    | imperfect   |  complex  | perfect   |imperfect  | complex
len(match) |                                      |
     1     |    *             *                   |
     2     |                  *              *    |      *           *
     3+    |                                 *    |                             *
        
     
         '''
        '''  
        #alignments_per_match
        #if closer to one, more likely to be global match
        num_alignments = 0
        for m in self.matches:
            num_alignments += len(m)
        alignments_per_match = num_alignments/len(self.matches)
        '''
        
        #selects a representative contig (the one with most matches) to test match alignment orientations
        labels = []
        
        repseq = self.get_rep_seq()
        
        #TODO - finding top_n has been moved into rep_seq function. Make sure that is connected in here to pass into label_finder
        query_matches = [m for m in self.matches if repseq in m.seqs]
        
        for m in query_matches:
            labels.append(m.label())
        
        #TODO - continue this function
        #uncovered regions 
        #can only be done within each match
        
        #end coverage? hold one contig and check alignent stats of all others against it?
        #would be expensive to do but would give confidence as to where alignments tend to happen.
        
        
        return labels
    
    def has_larger(self, graph):
        return
    
    def find_larger(self, graph):
        #TODO - not functioning yet
        '''find clusters that might envelop the sequences in the given cluster'''
        larger_components = graph.identify_larger(self)
        
        return larger_components
        
    def find_fragments(self, frag_matches):
        #TODO - is this in use?
        '''find fragments that constitute the same sequence in a truncated assembly'''
        frag_elements = []
        for m in frag_matches:
            if m.seqs[0] in self.nodes or m.seqs[1] in self.nodes:
                frag_elements.append(m)
        
        return frag_elements


'''---------------------End class definition-------------------------------'''

def cluster_agglomerate(cluster_objs, fragment_matches):
    '''connect clusters with other similar ones. connect with fragments'''
    #TODO - pipe for cluster_graph goes here
    return

def summary_file(cluster_objs, outfile):
    with open(outfile+".csv", mode="w") as f:
        writer = csv.writer(f, delimiter=",", quoting = csv.QUOTE_NONE)
        writer.writerow(['Size', 'Length', 'Coverage'])
        for c in cluster_objs:
            writer.writerow([c.size, c.av_length, c.av_cov])
        

def sort_clusters(cluster_list): #or maybe a dictionary instead?
    '''sort clusters by: 
    1. N in cluster, 2. length of N, 3. coverage of N
    Don't get it just from the names. That's SPAdes format but some people won't use spades
    Get it from the data directly. But megahit doesn't include any of this informaiton in the contig name, 
    and coverage would need reads mapped'''
    sortbysize = lambda c: (c.size is not None, c.size)
    sortbylength = lambda c: (c.av_length is not None, c.av_length)
    sortbycov = lambda c: (c.av_cov is not None, c.av_cov)
    
    
    #main sort function
    cluster_list.sort(reverse=True, key= lambda c: (sortbysize(c), sortbylength(c), sortbycov(c)))

def build_sig_match_dict(sig_matches): #to build dictionary of nodes linked to their match objects
    sig_match_dict = {}
    for m in sig_matches:
        
        link = m.seqs
        try:
            sig_match_dict[link[0]].append(m)
        except KeyError:
            sig_match_dict[link[0]] = [m]
        try:
            sig_match_dict[link[1]].append(m)
        except KeyError:
            sig_match_dict[link[1]] = [m]
            
    return sig_match_dict
    
def cluster_nucmer_matches(sig_matches): #sig_matches is a list of Nucmer_Match objects

    cluster_objs = []

    #TODO - maybe change this so that it only creates a sig_match dict for the most represented node in the cluster? Would save a lot of time in the sig_match_dict generation step.
    #which may be why the script is taking so long (>1 day currently)    
    sig_match_dict = build_sig_match_dict(sig_matches)
    
    cluster_list = union_find_cluster.union_find_pipe(sig_matches) #main clustering step
    
    for cluster in cluster_list:
        if len(cluster) >2:
            for node in cluster:
                nucmer_match_in_cluster = set(sig_match_dict[node])
            cluster_objs.append(Contig_Cluster(cluster, nucmer_match_in_cluster))

    return cluster_objs