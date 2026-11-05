#!/usr/bin/env bash
# prefetch + fasterq-dump + pigz
# Run from the terminal, in order. Variables P, IDX, T come from 02_setup.sh
cd $P/fastq
for run in $(tail -n +2 ../rna_samples.csv | cut -d, -f1); do
prefetch $run
fasterq-dump $run/$run.sra --split-files -e 8 -O .
pigz -p 8 ${run}_1.fastq ${run}_2.fastq
done
cd $P
ls -lh fastq/*.fastq.gz
