echo "GSE11431_Esrrb
GSE11431_Klf4
GSE11431_Oct4
GSE37262_Ncoa3
GSE11431_GFP
GSE11431_Nanog
GSE11431_Sox2" > samfiles

#!/bin/bash
while read LINE
do
/data/seqtools/samtools-1.1/samtools view -bS -o $LINE.bam $LINE.sam
done < /data/emmabell42/Coverage_plots/BAM_files/samfiles