#!/bin/bash


#remember that the assembly name is delimited from the node by __ (double underscore) this will need to go into repeatM
module load parallel

function amend_names {
assembly=$(basename $1)
sed  's/>/&'"$assembly"'__/' $1 > assembly_amended.$1
echo $1 amended
}

ls $@ > source_assemblies.txt

parallel -a source_assemblies.txt amend_names {}

#previous method
#for file in $@; do assembly=$(basename $file); sed 's/>/&'"$assembly"'__/' $file >> nucmer_ready.fa; echo $file "amended" ; done

wait

cat *assembly_amended* >> nucmer_ready.fa
