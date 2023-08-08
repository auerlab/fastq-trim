#!/bin/sh -e

fqt=SRR1972918_1-trimmed-smart10.fastq.zst
cutadapt=SRR1972918_1-trimmed-cutadapt.fastq.zst
trimmomatic=SRR1972918_1-trimmed-trimmomatic.fastq.zst

for file in $fqt $cutadapt $trimmomatic; do
    zstdcat $file | awk 'NR % 4 == 2' > $file-seqs.txt
done
wc -l *-seqs.txt

echo CTGTCTCTTATA
diff $fqt-seqs.txt $cutadapt-seqs.txt | more
echo CTGTCTCTTATA

printf "\nFastQ-Trim vs cutadapt:\n\n"
printf "Sequences in both files that differ:  "
diff -y $fqt-seqs.txt $cutadapt-seqs.txt | fgrep '|' | wc -l
printf "Sequences only in fastq-trim output:  "
diff -y $fqt-seqs.txt $cutadapt-seqs.txt | fgrep '<' | wc -l
printf "Sequences only in cutadapt output:    "
diff -y $fqt-seqs.txt $cutadapt-seqs.txt | fgrep '>' | wc -l

printf "\nTrimmomatic vs cutadapt:\n\n"
printf "Sequences in both files that differ:  "
diff -y $trimmomatic-seqs.txt $cutadapt-seqs.txt | fgrep '|' | wc -l
printf "Sequences only in trimmomatic output: "
diff -y $trimmomatic-seqs.txt $cutadapt-seqs.txt | fgrep '<' | wc -l
printf "Sequences only in cutadapt output:    "
diff -y $trimmomatic-seqs.txt $cutadapt-seqs.txt | fgrep '>' | wc -l

printf "\nFastQ-Trim vs Trimmomatic:\n\n"
printf "Sequences in both files that differ:  "
diff -y $fqt-seqs.txt $trimmomatic-seqs.txt | fgrep '|' | wc -l
printf "Sequences only in fastq-trim output:  "
diff -y $fqt-seqs.txt $trimmomatic-seqs.txt | fgrep '<' | wc -l
printf "Sequences only in trimmomatic output: "
diff -y $fqt-seqs.txt $trimmomatic-seqs.txt | fgrep '>' | wc -l

