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
./cave-man-install.sh
cd Test

export GZIP=-1

. ./fetch-infiles.sh

# After fetch-files.sh
case $1 in
big)
    suffix=''
    infile1=$long1
    infile2=$long2
    ;;

little)
    suffix='-250k'
    infile1=$short1
    infile2=$short2
    ;;

*)
    usage
    ;;

esac
shift

# adapter=AGATCGGAAGAGCACAC # Our mouse data
adapter=CTGTCTCTTATA
part=CTGTCTCTT
rand=TCGAACGGC
strict_percent=20


# Use gzip -1 for output to avoid bottleneck
outfile1_raw=${sample}_1$suffix-trimmed-raw.fastq.gz
outfile2_raw=${sample}_2$suffix-trimmed-raw.fastq.gz
outfile1_exact=${sample}_1$suffix-trimmed-exact.fastq.gz
outfile1_smart10=${sample}_1$suffix-trimmed-smart10.fastq.gz
outfile2_smart10=${sample}_2$suffix-trimmed-smart10.fastq.gz
outfile1_smart_strict=${sample}_1$suffix-trimmed-smart$strict_percent.fastq.gz
outfile1_cutadapt=${sample}_1$suffix-trimmed-cutadapt.fastq.gz
outfile2_cutadapt=${sample}_2$suffix-trimmed-cutadapt.fastq.gz
outfile1_trimmo=${sample}_1$suffix-trimmed-trimmomatic.fastq.gz
outfile1_paired=${sample}_1$suffix-trimmed-paired.fastq.gz
outfile2_paired=${sample}_2$suffix-trimmed-paired.fastq.gz
outfile1_fastp=${sample}_1$suffix-trimmed-fastp.fastq.gz

printf "\ngzip flags = $GZIP\n"

# Make sure all runs benefit equally from read buffering
printf "\nBuffering inputs...\n"
cat $infile1 $infile2 > /dev/null

infile1_uc=${infile1%.xz}
outfile1_uc=${sample}_1-out.fastq
xzcat $infile1 > $infile1_uc

printf "\ncat $infile1_uc > $outfile1_uc\n"
time cat $infile1_uc > $outfile1_uc
rm -f $outfile1_uc

printf "\nTrimming with uncompressed input and output...\n"
time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    $infile1_uc $outfile1_uc

printf "\nCutadapt 2 core uncompressed input and output...\n"
time cutadapt --report=minimal \
   --cores=2 --quality-cutoff=20 --minimum-length=30 -a $adapter \
   -o $outfile1_cutadapt $infile1_uc 2>&1 | fgrep -v reads/min

printf "\nTrimming with compressed input and uncompressed output...\n"
time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    $infile1 $outfile1_uc

printf "\nAll remaining tests use compressed input and output...\n"

printf "\nTiming compressed read and write without trimming...\n"
time xzcat $infile1 | gzip > $outfile1_raw
time xzcat $infile2 | gzip > $outfile2_raw

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    --exact-match \
    $infile1 $outfile1_exact

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    $infile1 $outfile1_smart10

time ../fastq-trim "$@" \
    --3p-adapter1 $adapter \
    --max-mismatch-percent $strict_percent \
    $infile1 $outfile1_smart_strict

# More than 1 core doesn't seem to help
for cores in 1; do
    # FIXME: This is not yet verified for comparison to other tools
    printf "\nfastp $cores thread...\n"
    time xzcat $infile1 \
	| fastp -q 20 -l 30 --adapter_sequence CTGTCTCTTATA \
	    --thread $cores --stdin -o $outfile1_fastp
done

for cores in 1 2; do
    printf "\nCutadapt $cores core...\n"
    time cutadapt --report=minimal \
       --cores=$cores --quality-cutoff=20 --minimum-length=30 -a $adapter \
       -o $outfile1_cutadapt $infile1 2>&1 | fgrep -v reads/min
done

if which trimmomatic > /dev/null 2>&1; then
    # Last number after .fa file is simple clip threshold.  A perfect match of
    # 12 bases give a score of about 7 according to docs.  5 was chosen by trial
    # and error to bring the missed partial adapters to a level similar to
    # fastq-trim and cutadapt.
    printf "\n"
    args="SE /dev/stdin $outfile1_trimmo ILLUMINACLIP:nextera.fa:2:30:5 TRAILING:20 MINLEN:30"
    time sh -c "xzcat $infile1 | trimmomatic $args"
    
    time ../fastq-trim "$@" \
	--3p-adapter1 $adapter \
	$infile2 $outfile2_smart10
fi

printf "\nCutadapt 2 core reverse read...\n"
time cutadapt --report=minimal \
   --cores=2 --quality-cutoff=20 --minimum-length=30 -a $adapter \
   -o $outfile2_cutadapt $infile2 2>&1 | fgrep -v reads/min

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

printf "Smart match $strict_percent output %-12s: " $part
gzcat $outfile1_smart_strict | fgrep $part | wc -l
printf "Smart match $strict_percent output %-12s: " $rand
gzcat $outfile1_smart_strict | fgrep $rand | wc -l

printf "Cutadapt output %-12s:       " $part
gzcat $outfile1_cutadapt | fgrep $part | wc -l
printf "Cutadapt output %-12s:       " $rand
gzcat $outfile1_cutadapt | fgrep $rand | wc -l

if which trimmomatic > /dev/null 2>&1; then
    printf "Trimmomatic output %-12s:    " $part
    gzcat $outfile1_trimmo | fgrep $part | wc -l
    printf "Trimmomatic output %-12s:    " $rand
    gzcat $outfile1_trimmo | fgrep $rand | wc -l
fi

