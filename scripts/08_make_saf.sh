#!/usr/bin/env bash
# Gene annotation in SAF format
# Run from the terminal, in order. Variables P, IDX, T come from 02_setup.sh
Rscript -e 'suppressMessages(library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)); g <- genes(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene); write.table(data.frame(GeneID=g$gene_id, Chr=as.character(seqnames(g)), Start=start(g), End=end(g), Strand=as.character(strand(g))), "counts/genes.saf", sep="\t", quote=FALSE, row.names=FALSE)'
head counts/genes.saf
