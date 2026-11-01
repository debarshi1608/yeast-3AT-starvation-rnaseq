#!/bin/bash
set -e
THREADS=8
GENOME_INDEX="ref/sacCer3/sacCer3"
GTF="ref/sacCer3/Saccharomyces_cerevisiae.R64-1-1.105.gtf"
RAW_DIR="data/fastq"
ALIGN_DIR="results/aligned"
COUNTS_DIR="results/counts"

mkdir -p $ALIGN_DIR $COUNTS_DIR logs
echo "=== Starting RNA-Seq Preprocessing Pipeline ==="

for R1 in ${RAW_DIR}/*_1.fastq.gz; do
    BASE=$(basename "$R1" _1.fastq.gz)
    R2="${RAW_DIR}/${BASE}_2.fastq.gz"
    echo "Processing sample: $BASE"
    hisat2 -p $THREADS -x $GENOME_INDEX -1 "$R1" -2 "$R2" \
        2> "logs/${BASE}_hisat2.log" | \
        samtools view -bS - | \
        samtools sort -@ $THREADS -o "${ALIGN_DIR}/${BASE}.sorted.bam"
    samtools index "${ALIGN_DIR}/${BASE}.sorted.bam"
done

echo "=== Running featureCounts ==="
featureCounts -T $THREADS -p -a $GTF -o "${COUNTS_DIR}/gene_counts.txt" ${ALIGN_DIR}/*.sorted.bam
echo "=== Pipeline Completed Successfully ==="
