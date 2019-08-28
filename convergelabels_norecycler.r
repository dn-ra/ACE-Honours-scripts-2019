#import data.table module
library(data.table)
## strip arguments from command line
args = commandArgs(trailingOnly = TRUE)

##error if not 3 input files
if (length(args)!=3){
  stop('expected 3 input files: plasflow_predictions.tsv, cbar_results.txt, virsorter_output')
}

#load data from 3 methods
plasdata <- args[which(grepl('plasflow_predictions', args))]
cbardata <- args[which(grepl('cBar_result', args))]
virsortdata <- args[which(grepl('VIRSorter', args))]

##read data
print('reading plasflow')
plasflow <- read.table(plasdata,header=TRUE)
print('reading cbar')
cbar <- read.table(cbardata, sep="\t", col.names=c("contig_name", "contig_length", "label"))
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



#virsorter- needs to be VIRSorter_global-phage-signal.csv
virsort <- as.character(virsort$V1)
virsort <- sub(pattern = "VIRSorter_", replacement = "", x = virsort)
subnode <- substring(virsort, first = 1, last = 20)
posnodes <- pmatch(subnode,labeltable$contig_name)
labeltable$virsorter <- matrix(NA, nrow(plasflow))
labeltable$virsorter[posnodes] <- "phage"


##output
print('table construction complete!')
write.table(labeltable, file = 'predictions.tsv', sep ="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

#check if each of the three prediciton tools identify a single contig as plasmid
labeltable_dt <- data.table(labeltable)
plasmidlist <- labeltable_dt[plasflow_label=='plasmid' & cbar_label=='plasmid' & is.na(virsorter)]
cnt <- plasmidlist[, .N]

#export 3-way agreed plasmid prediction table
write.table(plasmidlist, file= 'predicted_plasmids.csv', row.names = FALSE, col.names = TRUE, quote=FALSE)

#final output to console
sprintf('%i contigs have three-way agreement on plasmid classification using cbar, plasflow, and virsorter', cnt)
print('table written to \'predicted_plasmids.csv\'')




