#!/bin/sh -e

make clean all

time ./fastq-trim --exact-match 250k-R1.fastq.xz 250k-R1-trimmed-ft.fastq

time ./fastq-trim 250k-R1.fastq.xz 250k-R1-trimmed-ft-smart.fastq

# Use parameters to match defaults for fastq-trim
# 3 or 4 cores make very little difference
printf "cutadapt...\n"
time cutadapt \
    --cores=2 --quality-cutoff=20 --minimum-length=30 -a AGATCGGAAGAG \
    -o 250k-R1-trimmed-ca.fastq 250k-R1.fastq.xz | fgrep adapters

printf "Taking first 25k reads...\n"
head -100000 250k-R1-trimmed-ca.fastq > ca-short.fastq
head -100000 250k-R1-trimmed-ft.fastq > ft-short.fastq
head -100000 250k-R1-trimmed-ft-smart.fastq > ft-smart-short.fastq

printf "Diffing cutadapt vs fastq-trim default...\n"
diff ca-short.fastq ft-short.fastq | grep '^>' | egrep -v '@|\+|FF' | wc

printf "Diffing cutadapt vs fastq-trim smart match...\n"
diff ca-short.fastq ft-smart-short.fastq | grep '^>' | egrep -v '@|\+|FF' | wc

