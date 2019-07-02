#!/bin/bash



awk 'BEGIN{RS=">"} {print ">"$0 >> "cluster1_seq"FNR-1".fa"}' $1
