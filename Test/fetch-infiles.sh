#!/bin/sh -e

if [ ! -e $infile1 ]; then
    printf "Downloading FASTQ files...\n"
    fastq-dump --maxSpotId 1000000 --split-files $sample
    xz ${infile1%.xz} &
    xz ${infile2%.xz}
    wait
fi

if [ ! -e $short1 ]; then
    printf "Compressing...\n"
    xzcat -3 $infile1 | head -1000000 | xz > $short1 &
    xzcat -3 $infile2 | head -1000000 | xz > $short2
    wait
fi

