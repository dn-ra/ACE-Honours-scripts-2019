#!/bin/bash

for file in $@; do assembly=$(basename $file); sed 's/>/&'"$assembly"'__/' $file > nucmer_ready.fa; done
