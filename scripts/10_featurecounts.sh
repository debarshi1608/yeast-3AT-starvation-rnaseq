#!/usr/bin/env bash
# Count reads for all samples (-s 2)
# Run from the terminal, in order. Variables P, IDX, T come from 02_setup.sh
BAMS=$(tail -n +2 rna_samples.csv | cut -d, -f2 | sed 's#^#align/#; s#$#.bam#' | tr '\n' ' ')
featureCounts -T $T -p --countReadPairs -s 2 -B -C -a counts/genes.saf -F SAF \
-o counts/counts_raw.txt $BAMS 2> counts/featureCounts.log
column -t counts/counts_raw.txt.summary | cut -c1-150
