#!/bin/bash

#can't remember what this was for??

#first argument is csv of clustered contigs.
#second argument is ls of fasta files where these contigs are located. wildcard [*] is accepted


CLUSTERS=$1
FASTASOURCE=${@} #still struggling with how to get 2nd until last elements from the argument array!
i=0
echo $CLUSTERS

echo $FASTASOURCE

echo ${@}

while IFS="," read -a cols; do
	i=$((i+1))
	echo ${cols[*]} | sed 's/ /\n/g' | sed 's/^/>/' > newlineseq
	awk 'BEGIN{RS=">";FS="\n"}NR==FNR{a[$1]++}NR>FNR{if ($1 in a && $0!="") printf ">%s",$0}' newlineseq $FASTASOURCE > cluster_$i.fa 
	grep '>' cluster_$i.fa | wc -l
done < "$CLUSTERS"

echo processed $i clusters
