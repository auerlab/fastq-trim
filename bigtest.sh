#!/bin/sh -e

#############################################################################
# Timing an early version with just adapter removal:
#
# Input: chondro-sample1-rep1-time1-R1.fastq.xz    1.0G
#
# Reads: 32455547  Adapters: 983636  Qual < 20: 1596887  Len < 30: 159
# xz -2:    642.24 real       845.04 user         9.95 sys  1.5G
# xz -1:    470.41 real       671.50 user         9.16 sys  1.5G
# gzip -5:  562.98 real       761.33 user         7.57 sys  1.5G
# gzip -1:  128.48 real       331.70 user         9.30 sys
# None:     96.38 real       145.31 user         6.80 sys   8.3G
#
# gzip -1 provides a good load balance for xzipped input:
# 95.92% xzcat
# 95.40% gzip
# 82.94% fastq-tr

make clean all
export GZIP=-1

if [ ! -e chondro-sample1-rep1-time1-R1-trimmed.fastq.gz ]; then
    time ./fastq-trim "$@" --3p-adapter1 AGATCGGAAGAGCACAC \
	chondro-sample1-rep1-time1-R1.fastq.xz \
	chondro-sample1-rep1-time1-R1-trimmed.fastq.gz
fi

if [ ! -e chondro-sample1-rep1-time1-R1-trimmed-smart.fastq.gz ]; then
    time ./fastq-trim "$@" --3p-adapter1 AGATCGGAAGAGCACAC \
	--adapter-smart-match \
	chondro-sample1-rep1-time1-R1.fastq.xz \
	chondro-sample1-rep1-time1-R1-trimmed-smart.fastq.gz
fi

cat << EOM

Scanning for a portion of the adapter to flag any adapters missed due to
base substitutions.

Also scanning for random sequence of the same length for comparison.

EOM

set +e  # Don't terminate when fgrep finds no matches
printf "Raw data AGATCGGAAGAG:              "
xzcat chondro-sample1-rep1-time1-R1.fastq.xz | fgrep AGATCGGAAGAG | wc -l

printf "Raw data random:                    "
xzcat chondro-sample1-rep1-time1-R1.fastq.xz | fgrep CCTGAGATCTTC | wc -l

printf "Exact match output AGATCGGAAGAG:    "
gzcat chondro-sample1-rep1-time1-R1-trimmed.fastq.gz | fgrep AGATCGGAAGAG | wc -l

printf "Exact match output random:          "
gzcat chondro-sample1-rep1-time1-R1-trimmed.fastq.gz | fgrep CCTGAGATCTTC | wc -l

printf "Smart match output AGATCGGAAGAG:    "
gzcat chondro-sample1-rep1-time1-R1-trimmed-smart.fastq.gz | fgrep AGATCGGAAGAG | wc -l

printf "Smart match output random:          "
gzcat chondro-sample1-rep1-time1-R1-trimmed-smart.fastq.gz | fgrep CCTGAGATCTTC | wc -l
