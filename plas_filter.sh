#!/bin/bash

TAG_FILE=$1

echo "finding all plasmid labels"

awk -v SEP=$2  'BEGIN{RS=SEP; FS="\n"} 
{
if($0~/plasmid/)
print SEP $0}' $TAG_FILE

