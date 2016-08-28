#!/bin/sh

path="/Volumes/Phisan Segate/case-control-simulate/"
name=""
for i in {10..12}
do
    if [[ $i -lt 10 ]]; then
        name="000$i.csv"
    elif [[ $i -lt 100 ]]; then
        name="00$i.csv"
    else
        name="$i.csv"
    fi

    arg=$path$name
    Rscript skat.R "$arg"
done
