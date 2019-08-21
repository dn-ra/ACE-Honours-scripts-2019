#!/bin/bash
#first input is clustered ORF_labels file. Second input is a string used to delimit separate contig records


TAG_FILE=$1

if [ -z "$2" ]
	then
		echo "you must also provide a string argument that matches the beginning of the contig names (eg. \'NODE\' or \'k_141\')"
	exit 1
fi

echo "finding all plasmid labels"

awk -v SEP=$2  'BEGIN{RS=SEP; FS="\n"} 
{
if($0~/plasmid/)
print SEP $0}' $TAG_FILE

