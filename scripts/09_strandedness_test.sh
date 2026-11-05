#!/usr/bin/env bash
# featureCounts with -s 0 / 1 / 2
# Run from the terminal, in order. Variables P, IDX, T come from 02_setup.sh
for s in 0 1 2; do
featureCounts -T $T -p --countReadPairs -s $s -B -C -a counts/genes.saf -F SAF \
-o counts/strand_s$s.txt align/t00_r1.bam > /dev/null 2>&1
echo "-s $s assigned pairs: $(awk '$1=="Assigned"{print $2}' counts/strand_s$s.txt.summary)"
done
