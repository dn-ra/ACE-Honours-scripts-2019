#!/bin/bash

module -s load prokka
module -s load diamond


DELTAMATCHES=$(readlink -f $1)
FASTADIR=$(readlink -f $2)

echo $FASTADIR

filename=$(basename -- $DELTAMATCHES)
base=${filename%.*}
TMPDIR="match2prokka_"$base
mkdir -v "$TMPDIR"
cd $TMPDIR

#split matches into separate files, retain order
NM=2
echo processing $NM columns of matches

#exract contig names and set source fasta file varables
FASTAFILES=$(awk -v NM=$NM 'BEGIN{RS="\n";FS=" "} NR==1{print $1, $2}NR>1{for(i=1;i<=NM;i++){name="names_"i;print ">"$i> name}}' $DELTAMATCHES)
echo $FASTAFILES

#retrieve sequences from original fasta file

#for file in names_*; do
for ((i=1;i<=2;i++)); do

file="names_"$i
fasta=$(echo $FASTAFILES | cut -d " " -f $i)

echo names file is $file
echo source fastafile is $fasta

awk 'BEGIN{RS=">";FS="\n"}NR==FNR{a[$1]++}NR>FNR{if ($1 in a && $0!="") printf ">%s",$0}' $file  $FASTADIR/$fasta > $file"_seqs"

done


#prokka
#TODO - add in CARD and BACMET databases

for file in names_*_seqs; do


echo running $file through PROKKA
prokka $file --prefix "PROKKA_"$file --quiet

done

#diamond


echo cleaning up...
module unload prokka
module unload diamond


cd ..
