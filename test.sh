#!/bin/sh -e

make clean all
export GZIP=-1  # Speed up output compression

time ./fastq-trim "$@" \
    --3p-adapter1 AGATCGGAAGAGCACAC \
    250k-R1.fastq.xz 250k-R1-trimmed-ft.fastq.gz

time ./fastq-trim "$@" \
    --3p-adapter1 AGATCGGAAGAGCGTCG \
    250k-R2.fastq.xz 250k-R2-trimmed-ft.fastq.gz

printf "Smart adapter matching...\n"
time ./fastq-trim "$@" \
    --adapter-smart-match \
    --3p-adapter1 AGATCGGAAGAGCACAC \
    250k-R1.fastq.xz 250k-R1-smart-trimmed-ft.fastq.gz

printf "Paired mode...\n"
time ./fastq-trim "$@" \
    --3p-adapter1 AGATCGGAAGAGC \
    --3p-adapter2 AGATCGGAAGAGC \
    250k-R1.fastq.xz 250k-R1-paired-trimmed-ft.fastq.gz \
    250k-R2.fastq.xz 250k-R2-paired-trimmed-ft.fastq.gz

cat << EOM

Scanning for a portion of the adapter to flag any adapters missed due to
base substitutions.

Also scanning for random sequence of the same length for comparison.

EOM

set +e  # Don't terminate when fgrep finds no matches
printf "Raw data AGATCGGAAGAG:              "
xzcat 250k-R1.fastq.xz | fgrep AGATCGGAAGAG | wc -l
printf "Raw data random:                    "
gzcat 250k-R1-trimmed-ft.fastq.gz | fgrep CCTGAGATCTTC | wc -l

printf "Exact match output AGATCGGAAGAG:    "
gzcat 250k-R1-trimmed-ft.fastq.gz | fgrep AGATCGGAAGAG | wc -l

printf "Exact match output random:          "
gzcat 250k-R1-trimmed-ft.fastq.gz | fgrep CCTGAGATCTTC | wc -l

printf "Smart match output AGATCGGAAGAG:    "
gzcat 250k-R1-smart-trimmed-ft.fastq.gz | fgrep AGATCGGAAGAG | wc -l
printf "Smart match output random:          "
gzcat 250k-R1-trimmed-ft.fastq.gz | fgrep CCTGAGATCTTC | wc -l
