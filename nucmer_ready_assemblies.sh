#!/bin/bash
#edits contig names to include the assembly that they are sourced from
#for use in nucmer for repeatm

#remember that the assembly name is delimited from the node by __ (double underscore) this will need to go into repeatM
module load parallel

#tee this so that it sends to amended source file and to concatenated nucmer_ready.fasta

function amend_names {
assembly=$(basename $1)
echo $1 $assembly
sed  's/>/&'"$assembly"'__/' $1 > assembly_amended."$assembly" #stored in current directory
echo $1 amended
}
export -f amend_names


ls $@ > source_assemblies.txt



parallel --env amend_names -a source_assemblies.txt --verbose amend_names {}

#previous method
#for file in $@; do assembly=$(basename $file); sed 's/>/&'"$assembly"'__/' $file >> nucmer_ready.fa; echo $file "amended" ; done

wait

cat *assembly_amended* >> nucmer_ready.fa

module unload parallel
