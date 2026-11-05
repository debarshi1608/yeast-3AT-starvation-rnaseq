#!/usr/bin/env bash
# HISAT2 alignment, sorted BAM
# Run from the terminal, in order. Variables P, IDX, T come from 02_setup.sh
tail -n +2 rna_samples.csv | while IFS=, read run sample cond rep; do
hisat2 -p $T --max-intronlen 2000 -x $IDX \
-1 fastq/${run}_1.fastq.gz -2 fastq/${run}_2.fastq.gz \
--summary-file align/$sample.hisat2.txt 2>/dev/null \
| samtools sort -@ 4 -o align/$sample.bam -
samtools index align/$sample.bam
done
grep "overall alignment" align/*.hisat2.txt # alignment rates (Section 3)
samtools flagstat align/t00_r1.bam # optional extra check
