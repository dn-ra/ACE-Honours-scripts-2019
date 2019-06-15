#!/bin/bash

module -s load prokka
module -s load diamond

DELTAMATCHES=$(readlink -f $1)
FASTADIR=$(readlink -f $2)



#FASTA1=$(readlink -f $2)
#FASTA2=$(readlink -f $2)


filename=$(basename -- $DELTAMATCHES)
base=${filename%.*}
TMPDIR="match2prokka_"$base
mkdir -v "$TMPDIR"
cd $TMPDIR

#split matches into separate files, retain order
NM=$(awk 'NR==1{print NF}'  $DELTAMATCHES)
echo processing $NM columns of matches

awk -v NM=$NM 'BEGIN{RS="\n";FS=" "} {for(i=1;i<=NM;i++){name="names_"i;print ">"$i> name}}' $DELTAMATCHES


#retrieve sequences from original fasta file

for file in names_*; do

#FASTAFILE=


echo tmp file is $file

awk 'BEGIN{RS=">";FS="\n"}NR==FNR{a[$1]++}NR>FNR{if ($1 in a && $0!="") printf ">%s",$0}' $file  $FASTAFILE > $file"_seqs"

done


#prokka
#TODO - add in CARD and BACMET databases

for file in *_seqs; do


echo running $file through PROKKA
#prokka $file --prefix "PROKKA_"$file --quiet

done

#diamond


#elements stolen from
#https://stackoverflow.com/questions/24116307/split-file-into-multiple-files-by-columns

echo cleaning up...
module unload prokka
module unload diamond


cd ..
