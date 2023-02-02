#!/bin/sh -e

fqt=SRR1972918_1-trimmed-smart10.fastq.gz
cutadapt=SRR1972918_1-trimmed-cutadapt.fastq.gz

for file in $fqt $cutadapt; do
    zcat $file | awk 'NR % 4 == 2' > $file-seqs.txt
done
wc *-seqs.txt
printf "Sequences in both files that differ: "
diff -y *-seqs.txt | fgrep '|' | wc -l
printf "Sequences only in fastq-trim output: "
diff -y $fqt-seqs.txt $cutadapt-seqs.txt | fgrep '<' | wc -l
printf "Sequences only in cutadapt output:   "
diff -y $fqt-seqs.txt $cutadapt-seqs.txt | fgrep '>' | wc -l

