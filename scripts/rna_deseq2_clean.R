# 3-AT time-course RNA-seq, DESeq2 analysis (GSE110413). Run after 10_featurecounts.sh.
# Blocks follow the order of the report: steps 2.11-2.16 must run in order, then the plots.

# ---- DESeq2 step 1: libraries and settings ----
suppressMessages({library(DESeq2); library(ggplot2); library(pheatmap); library(ggrepel)})
setwd('~/3AT_RNAseq') # project folder that contains counts/ and rna_samples.csv
padj_cut <- 0.05 # adjusted p-value cut-off
lfc_cut <- 1 # |log2 fold change| >= 1 (2-fold)
dir.create('results', showWarnings = FALSE); dir.create('plots', showWarnings = FALSE)
cols <- c(up = '#c0392b', down = '#2471a3', ns = 'grey70') # colours used in the plots

# ---- step 2: read counts and sample sheet ----
cnt <- read.delim('counts/counts_raw.txt', comment.char = '#', check.names = FALSE)
mat <- as.matrix(cnt[, 7:ncol(cnt)]) # columns 1-6 are annotation, 7+ are samples
rownames(mat) <- cnt$Geneid
colnames(mat) <- sub('\\.bam$', '', basename(colnames(mat)))
meta <- read.csv('rna_samples.csv', stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
mat <- mat[, meta$sample] # same sample order in counts and metadata
lv <- c('t00', 't01', 't05', 't10', 't20', 't40')
meta$condition <- factor(meta$condition, levels = lv) # t00 = reference level
lab <- c(t01 = '1 min', t05 = '5 min', t10 = '10 min', t20 = '20 min', t40 = '40 min')

# ---- step 3: gene symbols ----
sym <- rownames(mat)
if (requireNamespace('org.Sc.sgd.db', quietly = TRUE)) {
s <- AnnotationDbi::mapIds(org.Sc.sgd.db::org.Sc.sgd.db, keys = rownames(mat),
column = 'COMMON', keytype = 'ORF', multiVals = 'first')
sym <- ifelse(is.na(s), rownames(mat), s)
}
names(sym) <- rownames(mat)

# ---- step 4: build object, filter, fit ----
dds <- DESeqDataSetFromMatrix(countData = mat, colData = meta, design = ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ] # keep genes with >=10 reads in >=2 samples
dds <- DESeq(dds) # size factors, dispersions, Wald test
vsd <- vst(dds, blind = FALSE) # variance-stabilised values for plots
cat('genes tested:', nrow(dds), '\n')
resultsNames(dds)

# ---- step 5: results for every time point ----
tp <- names(lab)
res <- lapply(tp, function(t) {
r <- results(dds, contrast = c('condition', t, 't00'), alpha = padj_cut)
sh <- lfcShrink(dds, coef = paste0('condition_', t, '_vs_t00'), type = 'apeglm', quiet = TRUE)
d <- data.frame(gene = rownames(r), symbol = sym[rownames(r)], baseMean = r$baseMean,
log2FC = r$log2FoldChange, log2FC_shrunk = sh$log2FoldChange,
pvalue = r$pvalue, padj = r$padj, stringsAsFactors = FALSE)
d$direction <- 'ns'
sig <- !is.na(d$padj) & d$padj < padj_cut & abs(d$log2FC_shrunk) >= lfc_cut
d$direction[sig & d$log2FC_shrunk > 0] <- 'up'
d$direction[sig & d$log2FC_shrunk < 0] <- 'down'
d$time <- lab[[t]]
d
})
names(res) <- tp

# ---- step 6: save result tables ----
for (t in tp) write.csv(res[[t]][order(res[[t]]$padj), ],
sprintf('results/DE_%s_vs_untreated.csv', t), row.names = FALSE)
A <- do.call(rbind, res); A$time <- factor(A$time, levels = lab)
summ <- as.data.frame.matrix(table(A$time, factor(A$direction, levels = c('up', 'down', 'ns'))))
summ <- cbind(comparison = paste(rownames(summ), 'vs untreated'), summ)
write.csv(summ, 'results/DEG_summary_counts.csv', row.names = FALSE); print(summ)
wide <- res[[1]][, c('gene', 'symbol', 'baseMean')]
for (t in tp) {
wide[[paste0('lfc_', t)]] <- res[[t]]$log2FC_shrunk
wide[[paste0('padj_', t)]] <- res[[t]]$padj
wide[[paste0('dir_', t)]] <- res[[t]]$direction
}
wide$n_timepoints_DE <- rowSums(wide[, paste0('dir_', tp)] != 'ns')
write.csv(wide[order(-wide$n_timepoints_DE), ], 'results/DEG_all_timepoints.csv', row.names = FALSE)
write.csv(counts(dds, normalized = TRUE), 'results/normalized_counts.csv')

# ---- 4.1 DEG counts per time point (Fig 1) ----
B <- as.data.frame(table(time = A$time, direction = factor(A$direction, levels = c('up', 'down', 'ns'))))
B <- B[B$direction != 'ns', ]
B$y <- ifelse(B$direction == 'down', -B$Freq, B$Freq)
p5 <- ggplot(B, aes(time, y, fill = direction)) + geom_col() +
geom_text(aes(label = Freq), vjust = ifelse(B$direction == 'down', 1.2, -0.3), size = 3.5) +
scale_fill_manual(values = cols) + theme_bw() +
labs(x = 'time in 3-AT', y = 'number of DEGs (down < 0 < up)',
title = paste0('DEGs vs untreated (padj<', padj_cut, ', |log2FC|>=', lfc_cut, ')'))
ggsave('plots/05_DEG_counts.png', p5, width = 7, height = 5, dpi = 300)

# ---- 4.2 PCA (Fig 2) ----
pc <- plotPCA(vsd, intgroup = 'condition', returnData = TRUE)
pv <- round(100 * attr(pc, 'percentVar'))
p1 <- ggplot(pc, aes(PC1, PC2, colour = condition, label = name)) + geom_point(size = 4) +
geom_text_repel(show.legend = FALSE) + theme_bw() +
xlab(paste0('PC1: ', pv[1], '% variance')) + ylab(paste0('PC2: ', pv[2], '% variance')) +
ggtitle('PCA of variance-stabilised counts')
ggsave('plots/01_PCA.png', p1, width = 7, height = 5.5, dpi = 300)

# ---- 4.3 Sample distances (Fig 3) ----
ann <- data.frame(condition = meta$condition, row.names = meta$sample)
pheatmap(as.matrix(dist(t(assay(vsd)))), annotation_col = ann, main = 'Sample distances',
filename = 'plots/02_sample_distance.png', width = 7, height = 6)

# ---- 4.4 MA plots (Fig 4) ----
p3 <- ggplot(A, aes(baseMean, log2FC_shrunk, colour = direction)) +
geom_point(size = 0.6, alpha = 0.6) + scale_x_log10() + scale_colour_manual(values = cols) +
geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = 2) + facet_wrap(~ time, nrow = 2) +
theme_bw() + labs(x = 'mean normalised count', y = 'shrunken log2 fold change vs untreated',
title = 'MA plots, 3-AT vs untreated')
ggsave('plots/03_MA_plots.png', p3, width = 11, height = 6.5, dpi = 300)

# ---- 4.5 Volcano plots (Fig 5) ----
V <- A[!is.na(A$padj), ]; V$nlp <- -log10(V$padj)
sig <- V[V$direction != 'ns', ]
top <- do.call(rbind, lapply(split(sig, droplevels(sig$time)), function(x) head(x[order(x$padj), ], 8)))
p4 <- ggplot(V, aes(log2FC_shrunk, nlp, colour = direction)) + geom_point(size = 0.6, alpha = 0.6) +
scale_colour_manual(values = cols) + geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2) +
geom_hline(yintercept = -log10(padj_cut), linetype = 2) + facet_wrap(~ time, nrow = 2) + theme_bw() +
labs(x = 'shrunken log2 fold change vs untreated', y = '-log10 adjusted p-value',
title = 'Volcano plots, 3-AT vs untreated') +
geom_text_repel(data = top, aes(label = symbol), size = 2.5, max.overlaps = 20, show.legend = FALSE)
ggsave('plots/04_volcano.png', p4, width = 11, height = 6.5, dpi = 300)

# ---- 4.6 Heatmap of DEGs (Fig 6) ----
deg <- unique(unlist(lapply(res, function(d) d$gene[d$direction != 'ns'])))
ord <- order(meta$condition)
pheatmap(assay(vsd)[deg, ord], scale = 'row', cluster_cols = FALSE, show_rownames = FALSE,
annotation_col = ann, main = paste(length(deg), 'DEGs (union), row z-score'),
filename = 'plots/06_heatmap_DEGs.png', width = 7, height = 9)

# ---- 4.7 UpSet plots (Fig 7, 8) ----
library(UpSetR)
for (dr in c('up', 'down')) {
L <- lapply(res, function(d) d$gene[d$direction == dr]); names(L) <- lab
L <- L[lengths(L) > 0]
png(sprintf('plots/07_upset_%s.png', dr), width = 2000, height = 1400, res = 250)
print(upset(fromList(L), nsets = length(L), order.by = 'freq',
mainbar.y.label = paste('shared', dr, 'genes')))
dev.off()
}

# ---- 4.8 Known Gcn4 targets (Fig 9) ----
kn <- c(HIS4 = 'YCL030C', HIS3 = 'YOR202W', ARG1 = 'YOL058W', ARG3 = 'YJL088W',
ARG4 = 'YHR018C', ARO3 = 'YDR035W', LYS1 = 'YIR034C', LEU4 = 'YNL104C')
K <- data.frame(gene = names(kn), min = 0, lfc = 0)
for (t in tp) {
d <- res[[t]][match(kn, res[[t]]$gene), ]
K <- rbind(K, data.frame(gene = names(kn), min = as.numeric(sub(' min', '', lab[[t]])),
lfc = d$log2FC_shrunk))
}
p8 <- ggplot(K, aes(min, lfc, colour = gene)) + geom_line() + geom_point() + theme_bw() +
labs(x = 'minutes in 3-AT', y = 'log2 fold change vs untreated',
title = 'Known Gcn4 targets (should increase)')
ggsave('plots/08_known_Gcn4_targets.png', p8, width = 7, height = 5, dpi = 300)

# ---- 4.9 Dispersion estimates (Fig 10) ----
png('plots/09_dispersion.png', width = 1400, height = 1100, res = 200)
plotDispEsts(dds)
dev.off()

