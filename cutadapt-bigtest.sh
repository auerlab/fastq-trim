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

make clean all
export XZ_OPT=-1
export GZIP=-1
cutadapt --help

time cutadapt --cores=2 --quality-cutoff=20 -a AGATCGGAAGAG \
     -o chondro-sample1-rep1-time1-R1-trimmed.fastq.gz \
     chondro-sample1-rep1-time1-R1.fastq.xz
