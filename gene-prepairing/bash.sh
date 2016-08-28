#!/bin/bash
for var in {1..500}
do
  echo $var
  /usr/local/bin/Rscript /Users/phisan/ResearchCode/SNPsR/R/gene-prepairing/gen-snps.R
done
