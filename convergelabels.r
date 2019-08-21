#import data.table module
library(data.table)
## strip arguments from command line
args = commandArgs(trailingOnly = TRUE)

##error if not 3 input files
if (length(args)!=4){
  stop('expected 4 input files: plasflow_predictions.tsv, cbar_results.txt, recycler_output.XXX')
}

#load data from 4 methods
plasdata <- args[which(grepl('plasflow_predictions', args))]
cbardata <- args[which(grepl('cBar_result', args))]
recyclerdata <- args[which(grepl('blast', args))]
virsortdata <- args[which(grepl('VIRSorter', args))]

##recycler must have blast output in csv format (-outfmt 10)
print('reading plasflow')
plasflow <- read.table(plasdata,header=TRUE)
print('reading cbar')
cbar <- read.table(cbardata, sep="\t", col.names=c("contig_name", "contig_length", "label"))
print('reading recycler')
recycler <- read.csv(recyclerdata, header=FALSE)
print('reading virsorter')
virsort <- read.csv(virsortdata, header = FALSE, comment.char = "#")

#plasflow - remove unwanted columns, parse classification strings
labeltable <- plasflow[-c(1, 4, 6:33)]
labeltable$label <- gsub("\\..*$", "", plasflow$label)
labeltable$contig_name <- as.character(plasflow$contig_name)
colnames(labeltable)[3] <- 'plasflow_label'

#cbar - match label case, treat as string
cbar$label <- sapply(cbar$label, tolower)
cbar$label <- as.character(cbar$label)
labeltable$cbar_label <- cbar$label


#recycler - add column in labeltable, extract contig names and mark as plasmid, keep info on blast hit length - aren't the full contigs!
labeltable$recycler_label <- matrix(NA, nrow(plasflow))
labeltable$recyclerblastlength <- matrix(NA, nrow(plasflow))
recycler$V2 <- as.character(recycler$V2)

#virsorter- needs to be VIRSorter_global-phage-signal.csv
virsort <- as.character(virsort$V1)
virsort <- sub(pattern = "VIRSorter_", replacement = "", x = virsort)
subnode <- substring(virsort, first = 1, last = 20)
posnodes <- pmatch(subnode,labeltable$contig_name)
labeltable$virsorter <- matrix(NA, nrow(plasflow))
labeltable$virsorter[posnodes] <- "phage"

#include blast hit length from recycler
#! is it correct to use 'sum' in ln42? blast finds several different sections of the contig that lines up with the node
#! change the for loop for something else. This is copying the recycler DF on every iteration. Hence why it's slow!
posnodes <- match(recycler$V2, labeltable$contig_name)
labeltable$recycler_label[posnodes] <- 'plasmid'

for (name in recycler$V2){
  row = which(labeltable$contig_name==name)
  labeltable$recyclerblastlength[row] <- sum(recycler[which(recycler$V2 ==name), 4])
  
  }

##output
print('table construction complete!')
write.table(labeltable, file = 'predictions.tsv', sep ="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#check if each of the three prediciton tools identify a single contig as plasmid
labeltable_dt <- data.table(labeltable)
plasmidlist <- labeltable_dt[plasflow_label=='plasmid' & cbar_label=='plasmid' & recycler_label=='plasmid' & is.na(virsorter)]
cnt <- plasmidlist[, .N]

#export 4-way agreed plasmid prediction table
write.table(plasmidlist, file= 'predicted_plasmids.csv', row.names = FALSE, col.names = TRUE, quote=FALSE)

#final output to console
sprintf('%i contigs have three-way agreement on plasmid classification', cnt)
print('table written to \'predicted_plasmids.csv\'')




