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

cd ..
make clean all
cd Test

export GZIP=-1

. fetch-infiles.sh

# Use gzip -1 for output to avoid bottleneck
outfile1_exact=${sample}_1-trimmed-exact.fastq.gz
outfile1_smart10=${sample}_1-trimmed-smart10.fastq.gz
outfile1_smart20=${sample}_1-trimmed-smart20.fastq.gz
outfile1_cutadapt=${sample}_1-trimmed-cutadapt.fastq.gz

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    --exact-match \
    $infile1 $outfile1_exact

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    --max-mismatch-percent 10 \
    $infile1 $outfile1_smart10

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    --max-mismatch-percent 20 \
    $infile1 $outfile1_smart20

for cores in 1 2; do
    printf "\nCutadapt $cores core...\n"
    time cutadapt --report=minimal \
       --cores=$cores --quality-cutoff=20 --minimum-length=30 -a $adapter \
       -o $outfile1_cutadapt $infile1 2>&1 | fgrep -v reads/min
done

cat << EOM

Scanning for a portion of the adapter to flag any adapters missed due to
base substitutions.

Also scanning for random sequence $rand for comparison.

EOM

set +e  # Don't terminate when fgrep finds no matches

printf "Raw data %-12s:              " $part
xzcat $infile1 | fgrep $part | wc -l

printf "Raw data %-12s:              " $rand
xzcat $infile1 | fgrep $rand | wc -l

printf "Exact match output %-12s:    " $part
gzcat $outfile1_exact | fgrep $part | wc -l

printf "Exact match output %-12s:    " $rand
gzcat $outfile1_exact | fgrep $rand | wc -l

printf "Smart match 10 output %-12s: " $part
gzcat $outfile1_smart10 | fgrep $part | wc -l

printf "Smart match 10 output %-12s: " $rand
gzcat $outfile1_smart10 | fgrep $rand | wc -l

printf "Smart match 20 output %-12s: " $part
gzcat $outfile1_smart20 | fgrep $part | wc -l

printf "Smart match 20 output %-12s: " $rand
gzcat $outfile1_smart20 | fgrep $rand | wc -l

printf "Cutadapt output %-12s:       " $part
gzcat $outfile1_cutadapt | fgrep $part | wc -l

printf "Cutadapt output %-12s:       " $rand
gzcat $outfile1_cutadapt | fgrep $rand | wc -l

