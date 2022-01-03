#!/bin/sh -e

#############################################################################
# Timing an early version with just adapter removal:
#
# Input: chondro-sample1-rep1-time1-R1.fastq.xz    1.0G
#
# Total reads processed = 32455547
# xz -2:    642.24 real       845.04 user         9.95 sys  1.5G
# xz -1:    470.41 real       671.50 user         9.16 sys  1.5G
# gzip -5:  562.98 real       761.33 user         7.57 sys  1.5G
# gzip -1:  112.32 real       295.73 user         8.69 sys  1.9G
# None:     96.38 real       145.31 user         6.80 sys   8.3G
#
# gzip -1 provides a good load balance for xzipped input:
# 95.92% xzcat
# 95.40% gzip
# 82.94% fastq-tr

make clean all
export XZ_OPT=-1
export GZIP=-1
time ./fastq-trim --3p-adapter AGATCGGAAGAGCACAC \
    chondro-sample1-rep1-time1-R1.fastq.xz \
    chondro-sample1-rep1-time1-R1-trimmed.fastq.gz
