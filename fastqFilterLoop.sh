#!/bin/bash

while read LINE
do
  cp $LINE "file_to_filter.txt"
  R CMD BATCH illuminaQualFilter.R
done < fastQfiles
