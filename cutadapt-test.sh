#!/bin/sh -e

#############################################################################
# Input: chondro-sample1-rep1-time1-R1.fastq.xz    1.0G
#
# gzip -1:  604.41 real       845.64 user         7.37 sys
#
# gzip -1 provides a good load balance for xzipped input:
# 95.92% xzcat
# 95.40% gzip
# 82.94% fastq-tr

export GZIP=-1
output=250k-R1-trimmed-ca.fastq.gz

time cutadapt --quality-cutoff=20 --minimum-length=30 -a AGATCGGAAGAG \
     -o $output 250k-R1.fastq.xz
blt fastx-stats $output
