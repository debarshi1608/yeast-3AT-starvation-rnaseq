#!/usr/bin/env bash
# Download sacCer3 and build the HISAT2 index (once)
# Run from the terminal, in order. Variables P, IDX, T come from 02_setup.sh
cd $P/ref
wget https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz
hisat2-build -p $T sacCer3.fa sacCer3_hisat
ls $IDX.1.ht2
cd $P
