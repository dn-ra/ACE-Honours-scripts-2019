#!/srv/sw/miniconda3/bin/python

#import relevant modules
import pysam
import argparse

parser = argparse.ArgumentParser(description='Find plasmids within a bam file.')
parser.add_argument('-inputFilePath', '-input', '-i', action = 'store', help = 'the bam input file')
parser.add_argument('-outputFilePath', '-output', '-o', action = 'store', help = 'the output')
args = parser.parse_args()
dict = vars(args)
inputFile = str(dict.get("inputFilePath"))
outputFile = str(dict.get("outputFilePath"))


#Write a new file
newFile = open(outputFile, "w")

#opens bam file for reading
samfile = pysam.AlignmentFile(inputFile, "rb")

#create dictionary to store names and length of all contigs
plasmidList = []

#iterate through lines to find all contigs
contigCount = 0
header = "contigName" + "\t" + "numberOfMatches" + "\t" + "matchingReads" #+ "\t" + "forwardReads" + "\t" + "reverseReads"
newFile.write(header + "\n")
for line in samfile:
    contigIdentifier = line.reference_name    
    readListA1 = []
    readListA2 = []
    readListB1 = []    
    readListB2 = []
    contigLength = samfile.lengths[samfile.get_tid(contigIdentifier)]
    #create list of contigs that are likely to be plasmids
    #checks that contig is > 1000bp
    if contigLength > 1000:
        #find reads associated with first 500bp of contig
        first500Reads = samfile.fetch(str(contigIdentifier), 1, 500)
        for read in first500Reads:
            queryName = read.query_name
            #exclude those reads that are forward mapped
            if read.is_reverse:
                #include those that are first reads
                if read.is_read1:            
                    readListA1.append(queryName)
                #include those that are second reads
                if read.is_read2:
                    readListB1.append(queryName)
        #find reads associated with last 500bp of contig
        last500Reads = samfile.fetch(str(contigIdentifier), (contigLength-500), contigLength)
        for read in last500Reads:
            queryName = read.query_name
            #include those that are forward mapped
            if read.mate_is_reverse:
                #include those that are first reads
                if read.is_read1:
                    readListB2.append(queryName)
                #include those that are second reads
                if read.is_read2:
                    readListA2.append(queryName)
        #starts count for number of matches
        count_readListA = 0
        count_readListB = 0
        
        #creates list for matching identifiers
        matches_readListA = []
        matches_readListB = []
        #for each entry in the first read lists, find if it matches to a corresponding second read
        for firstRead in readListA1:
            for secondRead in readListA2:
                if firstRead == secondRead:
                    count_readListA =+ 1
                    matches_readListA.append(firstRead)
        for firstRead in readListB1:
            for secondRead in readListB2:
                if firstRead == secondRead:
                    matches_readListB.append(firstRead)
        match_readListB_length = len(matches_readListB)
        match_readListA_length = len(matches_readListA)
        #print contig number and associated information
        newFile.write(contigIdentifier + '\t' + str(match_readListA_length + match_readListB_length) + '\t' + str(matches_readListA + matches_readListB))# + '\t' + str(readListB2 + readListA2) + '\t' + str (readListA1 + readListB1))
        if contigIdentifier not in plasmidList:
            if count_readListA or count_readListB >0:
                plasmidList.append(contigIdentifier)
    contigCount += 1


newFile.close      
samfile.close
