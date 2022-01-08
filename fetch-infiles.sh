#!/bin/sh -e

sample=SRR11180057

infile1=${sample}_1.fastq.xz
infile2=${sample}_2.fastq.xz

if [ ! -e $infile1 ]; then
    fastq-dump --maxSpotId 1000000 --split-files $sample
    xz ${infile1%.xz} &
    xz ${infile2%.xz}
    wait
fi

short1=${sample}_1-250k.fastq.xz
short2=${sample}_2-250k.fastq.xz
if [ ! -e $short1 ]; then
    xzcat $infile1 | head -1000000 | xz > $short1 &
    xzcat $infile2 | head -1000000 | xz > $short2
    wait
fi

# adapter=AGATCGGAAGAGCACAC # Our mouse data
adapter=CTGTCTCTTATA
part=CTGTCTCTTAT
rand=CCCGGCTGGTA
