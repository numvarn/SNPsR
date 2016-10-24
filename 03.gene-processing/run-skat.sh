#!/bin/sh

path="/Volumes/Sirikanlaya/case-control-replicated/0.0/"
name=""
for i in {1..100}
do
     name="$i.csv"
     arg=$path$name
     Rscript skat.R "$arg"
done
