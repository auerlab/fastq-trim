#!/bin/sh -e

#############################################################################
# Timing an early version with just adapter removal:
#
# FIXME: Rig this to download some public data instead of the chondro sample
#
# gzip -1 provides a good load balance for xzipped input:
# 95.92% xzcat
# 95.40% gzip
# 82.94% fastq-tr

make clean all
export GZIP=-1

# sample=neuro-sample10-rep1-time1
sample=2M

adapter=AGATCGGAAGAGC
time ./fastq-trim "$@" \
    --3p-adapter1 $adapter \
    --exact-match \
    $sample-R1.fastq.xz \
    $sample-R1-trimmed-exact.fastq.gz

time ./fastq-trim "$@" \
    --3p-adapter1 $adapter \
    --max-mismatch-percent 10 \
    $sample-R1.fastq.xz \
    $sample-R1-trimmed-smart10.fastq.gz

time ./fastq-trim "$@" \
    --3p-adapter1 $adapter \
    --max-mismatch-percent 20 \
    $sample-R1.fastq.xz \
    $sample-R1-trimmed-smart20.fastq.gz

printf "Running cutadapt...\n"
time cutadapt --report=minimal \
    --cores=2 --quality-cutoff=20 --minimum-length=30 -a $adapter \
    -o $sample-R1-trimmed-cutadapt.fastq.gz $sample-R1.fastq.xz \
    2>&1 | fgrep -v reads/min

cat << EOM

Scanning for a portion of the adapter to flag any adapters missed due to
base substitutions.

Also scanning for random sequence of the same length for comparison.

EOM

set +e  # Don't terminate when fgrep finds no matches

adapterpart=AGATCGGAAG
random=CCTGAGATCT

printf "Raw data %-12s:              " $adapterpart
xzcat $sample-R1.fastq.xz | fgrep $adapterpart | wc -l

printf "Raw data random:                    "
xzcat $sample-R1.fastq.xz | fgrep $random | wc -l

printf "Exact match output %-12s:    " $adapterpart
gzcat $sample-R1-trimmed-exact.fastq.gz | fgrep $adapterpart | wc -l

printf "Exact match output random:          "
gzcat $sample-R1-trimmed-exact.fastq.gz | fgrep $random | wc -l

printf "Smart match 10 output %-12s: " $adapterpart
gzcat $sample-R1-trimmed-smart10.fastq.gz | fgrep $adapterpart | wc -l

printf "Smart match 10 output random:       "
gzcat $sample-R1-trimmed-smart10.fastq.gz | fgrep $random | wc -l

printf "Smart match 20 output %-12s: " $adapterpart
gzcat $sample-R1-trimmed-smart20.fastq.gz | fgrep $adapterpart | wc -l

printf "Smart match 20 output random:       "
gzcat $sample-R1-trimmed-smart20.fastq.gz | fgrep $random | wc -l

printf "Cutadapt output %-12s:       " $adapterpart
gzcat $sample-R1-trimmed-cutadapt.fastq.gz | fgrep $adapterpart | wc -l

printf "Cutadapt output random:             "
gzcat $sample-R1-trimmed-cutadapt.fastq.gz | fgrep $random | wc -l

