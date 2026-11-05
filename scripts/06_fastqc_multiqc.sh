#!/usr/bin/env bash
# Raw read QC
# Run from the terminal, in order. Variables P, IDX, T come from 02_setup.sh
fastqc -t 4 -o qc/fastqc fastq/*.fastq.gz
multiqc qc/fastqc -o qc/multiqc
