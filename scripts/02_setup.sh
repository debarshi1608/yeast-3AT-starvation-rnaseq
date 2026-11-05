#!/usr/bin/env bash
# Variables, activate env, folders
# Run from the terminal, in order. Variables P, IDX, T come from 02_setup.sh
export P=$HOME/3AT_RNAseq
export IDX=$P/ref/sacCer3_hisat
export T=12
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq_de
mkdir -p $P/{ref,fastq,qc/fastqc,qc/multiqc,align,counts,results,plots}
cd $P
for x in hisat2 samtools featureCounts prefetch fasterq-dump fastqc multiqc pigz Rscript; do
command -v $x >/dev/null || echo "tool not found: $x"
done
