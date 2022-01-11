#!/bin/sh -e

##########################################################################
#   Synopsis:
#       ./test.sh big|little
#
#   Description:
#       Benchmark fastq-trim against cutadapt, reporting run times and
#       sequences removed
#
#   Arguments:
#       $1:     big = process whole fastq, little = 250k reads
#       
#   History:
#   Date        Name        Modification
#   2022-01-11  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 big|little [extra fastq-trim flags]\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# -lt 1 ]; then
    usage
fi

# SRR11180057 looks like it's already been trimmed.
# Only 138 Nextera adapters found in almost a million reads and most not
# near the 3' end.  Probably natural sequences.
# SRR1553607 has only 203445 reads
sample=SRR1972918

case $1 in
big)
    suffix=''
    ;;

little)
    suffix='-250k'
    ;;

*)
    usage
    ;;

esac
shift
infile1=${sample}_1$suffix.fastq.xz
infile2=${sample}_2$suffix.fastq.xz

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

# adapter=AGATCGGAAGAGCACAC # Our mouse data
adapter=CTGTCTCTTATA
part=CTGTCTCTT
rand=TCGAACGGC

# Use gzip -1 for output to avoid bottleneck
outfile1_exact=${sample}_1$suffix-trimmed-exact.fastq.gz
outfile1_smart10=${sample}_1$suffix-trimmed-smart10.fastq.gz
outfile1_smart20=${sample}_1$suffix-trimmed-smart20.fastq.gz
outfile1_cutadapt=${sample}_1$suffix-trimmed-cutadapt.fastq.gz
outfile1_trimmo=${sample}_1$suffix-trimmed-trimmomatic.fastq.gz
outfile1_paired=${sample}_1$suffix-trimmed-paired.fastq.gz
outfile2_paired=${sample}_2$suffix-trimmed-paired.fastq.gz

# Make sure all runs benefit equally from read buffering
printf "Buffering input files...\n"
cat $infile1 > /dev/null
cat $infile2 > /dev/null

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    --exact-match \
    $infile1 $outfile1_exact

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
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

args="SE /dev/stdin $outfile1_trimmo ILLUMINACLIP:nextera.fa:2:30:10 TRAILING:20 MINLEN:30"
time sh -c "xzcat $infile1 | trimmomatic $args"

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    $infile1 $outfile1_paired $infile2 $outfile2_paired

cat << EOM

Scanning for a portion of the adapter to flag any adapters missed due to
base substitutions.

Also scanning for random sequence $rand as a control.

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

printf "Trimmomatic output %-12s:    " $part
gzcat $outfile1_trimmo | fgrep $part | wc -l
printf "Trimmomatic output %-12s:    " $rand
gzcat $outfile1_trimmo | fgrep $rand | wc -l

