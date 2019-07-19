#!/bin/bash

#nucmer alignment of all assemblies in input
#input assemblies must be node-assembly matched (created by "nucmer_ready_assemblies.sh")

#input: list of n assemblies to be processed
#ouput: n-1 .delta files to be processed by repeatM

module load mummer/4.0###

ASSEMBLIES = ${@}

echo Assemblies input:
echo $ASSEMBLIES

mkdir nucmer_out
touch assemblies_iter_built.fa

QUERY_SIZE=0

for ((i=0; i<$#-1; i++)); do
CURR_REF=${@[$i]};
CURR_QUE=${@[$(expr $i + 1)]};
ITER=$(expr $i + 1);

SEQ_COUNT=$(grep '>' $CURR_QUE | wc -l);
QUERY_SIZE=$(expr $QUERY_SIZE + $SEQ_COUNT);

echo adding "$SEQ_COUNT" to concatenated query file in iteration "$ITER";
cat $CURR_QUE >> assemblies_iter_built.fa;


echo processing "$QUERY_SIZE" sequences in query file against reference file "CURR_REF";

time nucmer "$CURR_REF" assemblies_iter_built.fa -p iter_"$ITER";
done

rm assemblies_iter_built.fa
mv iter_*.delta nucmer_out
