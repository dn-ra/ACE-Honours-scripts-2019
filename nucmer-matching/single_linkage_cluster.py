'''To be used here for testing of single_contig clustering from delta match files'''

'''
    imports
        '''
        
import os
import re
import subprocess
import tempfile
from logging import warning #also import exceptions?


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
            self.av_cov = cov_total / spades_seqs
            self.av_length = length_total / spades_seqs
        else:
            warning('Non-Spades contig names not supported yet')
    #necessary?
    def has_spades(self):
        #check if there is at least one spades contig in node_list
        spades = False
        for c in self.nodes:
            if spadespattern.match(c):
                spades = True
                break
        return spades
    
    #TODO - for running plasflow or cbar?
    def plas_label(self):
        return None
    
    #TODO - mkdir for each cluster
    def retrieve_seqs(self, assembly_dir, outdir = None, outfile = 'cluster_out.fa', return_node_assembly_dict =False):
        #all I need are nodes and location of source fasta files? use sequence library in repeatM
        #pop out assembly number from start of contig? use as input the source fastafiles?
        #Beware leaked processes
        #can get assembly_dir from delta output?
        if outdir == None:
            outdir = 'CLUSTER_size_{}_avlen_{}_avcov_{}'.format(self.size, self.av_length, self.av_cov)
            
        print(outdir)
        os.mkdir(outdir)
        node_assembly_dict = {}
        tmp_node_fnas = {}
        for n in self.nodes:
            nodesplit = n.split("__")
            node_assembly_dict[">"+nodesplit[1]] = nodesplit[0] #remainder (1st entry) into assembly list
            
        cluster_out = open(outdir + "/" + outfile, 'w')
        #parallelise eventually. That's why i've written it to dictionaries first
        for node, assembly in node_assembly_dict.items():
            #seq_tmp = tempfile.NamedTemporaryFile(delete=False) #not necessary. Don't need to pass file to samtools view
            #tmp_node_fnas[node] = seq_tmp.name #not necessary. Don't need to pass file to samtools view.
            f = open("/".join([assembly_dir, assembly]))
            for i in f:
                if i.find(node)!=-1:
                    #seq_tmp.write(i.encode())
                    cluster_out.write(i)
                    wholeseq =False #flag to tell me if I've taken the whole sequence yet
                    while wholeseq ==False:
                        line = next(f)
                        if line.startswith('>'):
                            wholeseq = True
                        else:
                            #seq_tmp.write(line.encode()) #not necessary
                            cluster_out.write(line)
            #seq_tmp.close() # not necessary
        cluster_out.close()
        
        if return_node_assembly_dict == True:
            return node_assembly_dict #return this to pass into gen_minibam. Not particularly neat
        
        #delete tempfiles, delete = False means user must handle their deletion
        #[os.removeitems(f) for f in tmp_node_fnas.values()] # note necessary
        
#        previously proposed method
        #awk (retrieve sequence) $contigID $assembly_file
#        nodefind_awk = '''-v name=$node 'BEGIN{RS=">";FS="\n"}NR>1{if ($1~/name/) print ">"$0}\''''
#        CMD = 'while read -r node <&3 && read -r assembly <&4; do awk {} $assembly; done 3< <(printf "{}") 4< <(printf "{}")'.format(nodefind_awk, nodestring, assemblystring)
#        subprocess.call(['bash', '-c', CMD])
            
#what if I look for orientation of reads mapped from other assemblies? Shouldn't that be consistent too?
#have to call from within retreieve_seqs for tempfiles to be accessible
    def gen_minibam(node_assembly_dict, bam_location, outdir = 'minibam_out'): #pass self into this?
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
        
    #TODO -
    def label_cluster(self):
        '''label cluster as linear or circular based on alignment evidence
        Evidence includes: 
        - Do regions align globally or differentially? #no of matches in Nucmer_Match
        - do endings overlap?
        - where do the reads go? when do I generate bam files? #look at minibam
        -are there any regions that don't match?
        '''
        
        #alignments_per_match
        #if closer to one, more likely to be global match
        num_alignments = 0
        for m in self.matches:
            num_alignments += len(m)
        alignments_per_match = num_alignments/len(self.matches)
        
        #overlapping regions 
        #can only be done within each match. cannot
        #assume that regions are directly comparable between match objects
        #set as boolean?
        overlap = 0
        seen_nodes = []
        
        for n in self.nodes:
            current = n
            seen_nodes.append(current)
            check_matches = [match for match in self.matches if current in match.seqs]
            for m in check_matches:
                ref = m.seqs.index(current) + 1
                starts = getattr(m, 'hitstarts_{}'.format(ref))
                stops = getattr(m, 'hitstops_{}'.format(ref))
        
        #uncovered regions 
        #can only be done within each match
        
        #end coverage? hold one contig and check alignent stats of all others against it?
        #would be expensive to do but would give confidence as to where alignments tend to happen.
        
        
        
        return None



'''---------------------End class definition-------------------------------'''


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

def cluster_nucmer_matches(sig_matches): #sig_matches is a list of Nucmer_Match objects
    sig_match_dict = {}
    match_links = []
    cluster_objs = []
    
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
        cluster_objs.append(Contig_Cluster(cluster, nucmer_match_in_cluster))

    return cluster_objs

'''This from RepeatM Clusterer module - author wwood'''
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
    
