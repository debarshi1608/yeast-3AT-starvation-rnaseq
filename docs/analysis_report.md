# Time-Course RNA-Seq Differential Expression Analysis of the 3-AT Starvation Response in Saccharomyces cerevisiae (GSE110413)

## 1. Aim and Background
3-Amino-1,2,4-triazole (3-AT) is a competitive inhibitor of imidazoleglycerol-phosphate dehydratase (His3p). It causes acute histidine starvation and activates the General Amino Acid Control (GAAC) pathway, whose master transcription factor is Gcn4. This study models the temporal dynamics of yeast exposed to 3-AT across 1, 5, 10, 20, and 40 minutes relative to untreated controls (0 min)[cite: 8].

## 2. Dataset and Summary Metrics
* **Sequencing:** 8 paired-end Illumina libraries from *S. cerevisiae* (sacCer3 reference genome)[cite: 8].
* **Alignment Rates:** Ranged from 89.41% (t01_r1) to 96.18% (t00_r2) using HISAT2[cite: 8].
* **Strandedness:** Confirmed reversely stranded library protocol (`-s 2` via `featureCounts`)[cite: 8].
* **Filtering:** Tested 6,358 genes after retaining those with $\ge 10$ reads in $\ge 2$ samples[cite: 8].

## 3. Key Findings & Temporal Waves
* **1–5 Min:** Minimal transcriptional response; pre-transcriptional lag phase (0 DEGs at 1 min, 6 DEGs at 5 min)[cite: 8].
* **10–20 Min:** Rapid induction of Gcn4 targets (73 DEGs at 10 min; 288 DEGs at 20 min)[cite: 8].
* **40 Min:** Large-scale transcriptomic remodelling with **1,408 DEGs** (660 up, 748 down)[cite: 8].
* **Biological Validation:** Classic Gcn4 targets (*HIS3*, *HIS4*, *ARG1*, *ARG3*, *ARG4*, *ARO3*, *LYS1*, *LEU4*) showed expected strong upregulation trajectories[cite: 8].
