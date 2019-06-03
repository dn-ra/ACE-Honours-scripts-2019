#import data.table module
library(data.table)
## strip arguments from command line
args = commandArgs(trailingOnly = TRUE)

##error if not 3 input files
if (length(args)!=3){
  stop('expected 3 input files: plasflow_predictions.tsv, cbar_results.txt, recycler_output.XXX')
}

#load data from 3 methods
plasdata <- args[which(grepl('plasflow', args))]
cbardata <- args[which(grepl('cbar', args))]
recyclerdata <- args[which(grepl('blast', args))]


##recycler must have blast output in csv format (-outfmt 10)
plasflow <- read.table(plasdata,header=TRUE)
cbar <- read.table(cbardata, sep="\t", col.names=c("contig_name", "contig_length", "label"))
recycler = read.csv(recyclerdata, header=FALSE)

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

#include blast hit length from recycler
#! is it correct to use 'sum' in ln42? blast finds several different sections of the contig that lines up with the node
#! change the for loop for something else. This is copying the recycler DF on every iteration. Hence why it's slow!
for (name in recycler$V2){
  row = which(labeltable$contig_name==name)
  labeltable$recycler_label[row] <- 'plasmid'
  labeltable$recyclerblastlength[row] <- sum(recycler[which(recycler$V2 ==name), 4])
  
}

##output
print('table construction complete!')
write.table(labeltable, file = 'predictions.tsv', sep ="\t", row.names=FALSE, col.names=TRUE)

#check if each of the three prediciton tools identify a single contig as plasmid
labeltable_dt <- data.table(labeltable)
plasmidlist <- labeltable_dt[plasflow_label=='plasmid' & cbar_label=='plasmid' & recycler_label=='plasmid']
cnt <- plasmidlist[, .N]

write.table(plasmidlist, file= '3-way_predicted_plasmids.csv', row.names = FALSE, col.names = TRUE, quote=FALSE)

sprintf('%i contigs have three-way agreement on plasmid classification', cnt)
print('table written to \'3 way_predicted_plasmids.csv\'')




