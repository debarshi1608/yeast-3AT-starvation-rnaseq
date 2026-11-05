#!/usr/bin/env bash
# Write rna_samples.csv
# Run from the terminal, in order. Variables P, IDX, T come from 02_setup.sh
cat > rna_samples.csv <<'EOF'
run,sample,condition,replicate
SRR6703656,t00_r1,t00,1
SRR6703662,t00_r2,t00,2
SRR6703657,t01_r1,t01,1
SRR6703658,t05_r1,t05,1
SRR6703659,t10_r1,t10,1
SRR6703660,t20_r1,t20,1
SRR6703661,t40_r1,t40,1
SRR6703663,t40_r2,t40,2
EOF
