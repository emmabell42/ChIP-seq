#!/bin/bash
ls -a -1 *.sam > samfiles
sed -i 's/.sam//g' samfiles

while read LINE
do
/data/seqtools/samtools-1.1/samtools view -bS -o $LINE.bam $LINE.sam
done < ./samfiles
