library(DESeq2)
library(tidyverse)
library(pheatmap)

counts_raw <- read.table("results/counts/gene_counts.txt", header=TRUE, row.names=1, sep="\t")
countData <- counts_raw[, 6:ncol(counts_raw)]
colnames(countData) <- gsub(".sorted.bam", "", colnames(countData))
colnames(countData) <- gsub(".*\\/", "", colnames(countData))

sampleNames <- colnames(countData)
colData <- data.frame(
  row.names = sampleNames,
  condition = factor(rep(c("Control", "Treat_3AT"), each=3)),
  time = factor(rep(c("0min", "30min", "60min"), times=2))
)
colData <- colData[colnames(countData), ]

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition + time + condition:time)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

res <- results(dds)
resOrdered <- res[order(res$pvalue), ]
write.csv(as.data.frame(resOrdered), file="results/deseq2_all_results.csv")

vsd <- vst(dds, blind=FALSE)
pdf("results/pca_plot.pdf", width=6, height=5)
plotPCA(vsd, intgroup=c("condition", "time"))
dev.off()
print("DESeq2 analysis completed. Outputs saved to results/")
