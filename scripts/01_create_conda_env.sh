#!/usr/bin/env bash
# Create the conda environment (once)
# Run from the terminal, in order. Variables P, IDX, T come from 02_setup.sh
source ~/miniforge3/etc/profile.d/conda.sh
conda create -y -n rnaseq_de -c conda-forge -c bioconda \
sra-tools hisat2 samtools subread fastqc multiqc pigz r-base \
bioconductor-deseq2 bioconductor-apeglm \
bioconductor-txdb.scerevisiae.ucsc.saccer3.sgdgene bioconductor-org.sc.sgd.db \
r-ggplot2 r-pheatmap r-ggrepel r-upsetr
