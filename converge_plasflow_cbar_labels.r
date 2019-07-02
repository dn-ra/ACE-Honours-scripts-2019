#import data.table module
library(data.table)
## strip arguments from command line
args = commandArgs(trailingOnly = TRUE)

##error if not 2 input files
if (length(args)!=2){
  stop('expected 2 input files: plasflow_predictions.tsv, cbar_results.txt')
}

#load data from 2 methods
plasdata <- args[which(grepl('plasflow', args))]
cbardata <- args[which(grepl('cbar', args))]


##recycler must have blast output in csv format (-outfmt 10)
plasflow <- read.table(plasdata,header=TRUE)
cbar <- read.table(cbardata, sep="\t", col.names=c("contig_name", "contig_length", "label"))

#plasflow - remove unwanted columns, parse classification strings
labeltable <- plasflow[-c(1, 4, 6:33)]
labeltable$label <- gsub("\\..*$", "", plasflow$label)
labeltable$contig_name <- as.character(plasflow$contig_name)
colnames(labeltable)[3] <- 'plasflow_label'

#cbar - match label case, treat as string
cbar$label <- sapply(cbar$label, tolower)
cbar$label <- as.character(cbar$label)
labeltable$cbar_label <- cbar$label



##output
print('table construction complete!')
write.table(labeltable, file = 'predictions.tsv', sep ="\t", row.names=FALSE, col.names=TRUE)

#check if each of the three prediciton tools identify a single contig as plasmid
labeltable_dt <- data.table(labeltable)
plasmidlist <- labeltable_dt[plasflow_label=='plasmid' & cbar_label=='plasmid']
cnt <- plasmidlist[, .N]

write.table(plasmidlist, file= '2-way_predicted_plasmids.csv', row.names = FALSE, col.names = TRUE, quote=FALSE)

sprintf('%i contigs have two-way agreement on plasmid classification', cnt)
print('table written to \'2 way_predicted_plasmids.csv\'')




