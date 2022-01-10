#!/bin/sh -e

cd ..
make clean all
cd Test

. fetch-infiles.sh

# Use gzip -1 for output to avoid bottleneck
outfile1_exact=${sample}_1-trimmed-exact.fastq.gz
outfile1_smart=${sample}_1-trimmed-smart.fastq.gz
outfile2_smart=${sample}_2-trimmed-smart.fastq.gz
outfile1_cutadapt=${sample}_1-trimmed-cutadapt.fastq.gz

export GZIP=-1  # Speed up output compression

time ../fastq-trim "$@" \
    --exact-match \
    --3p-adapter1 $adapter \
    $short1 $outfile1_exact

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    $short1 $outfile1_smart

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    --3p-adapter2 $adapter \
    $short1 $outfile1_smart $short2 $outfile2_smart

cat << EOM

Scanning for a portion of the adapter to flag any adapters missed due to
base substitutions.

Also scanning for random sequence of the same length for comparison.

EOM

set +e  # Don't terminate when fgrep finds no matches
printf "Raw data $part:              "
xzcat $short1 | fgrep $part | wc -l
printf "Raw data $rand:              "
xzcat $short1 | fgrep $rand | wc -l

printf "Exact match output $part:    "
zcat $outfile1_exact | fgrep $part | wc -l
printf "Exact match output $rand:    "
zcat $outfile1_exact | fgrep $rand | wc -l

printf "Smart match output $part:    "
zcat $outfile1_smart | fgrep $part | wc -l
printf "Smart match output $rand:    "
zcat $outfile1_smart | fgrep $rand | wc -l
