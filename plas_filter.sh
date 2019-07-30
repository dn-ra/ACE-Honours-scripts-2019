#!/bin/bash

TAG_FILE=$1

re='^[0-9]+$'
#if decimcal value included in $2, use it as threshold, else produce all contigs with at least 1 plasmid label

if [ -z $2 ]; then
echo "finding all plasmid labels"
THRESH=0;

elif (! [[ $2 =~ $re ]]) ; then
echo "error: Not a number"

else
echo "finding contigs with tag ratio threshold of $2 percent"
THRESH=$2
echo $THRESH
fi


awk -v THRESH=$THRESH 'BEGIN{RS="NODE"; FS="\n"} 
{VAR1="NODE" $1 "\n";PLASCOUNT=0.0;LABELCOUNT=0.0;

for (i=2; i<=NF; i++) 
{
if($i~/plasmid/) {VAR1=VAR1 $i "\n";PLASCOUNT+=1.0;LABELCOUNT+=1.0;} 
else if($i~/genomic/ || ($i~/chromosome/)) {LABELCOUNT+=1.0;}
} 
if (LABELCOUNT>0)
{ PERCENT=PLASCOUNT*100/LABELCOUNT
if ( PERCENT>THRESH )
print VAR1, PLASCOUNT, LABELCOUNT}}' $TAG_FILE

