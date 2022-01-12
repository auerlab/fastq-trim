#!/bin/sh -e

# SRR11180057 looks like it's already been trimmed.
# Only 138 Nextera adapters found in almost a million reads and most not
# near the 3' end.  Probably natural sequences.
# SRR1553607 has only 203445 reads
sample=SRR1972918

long1=${sample}_1$suffix.fastq.xz
long2=${sample}_2.fastq.xz

if [ ! -e $long1 ]; then
    printf "Downloading FASTQ files...\n"
    fastq-dump --maxSpotId 1000000 --split-files $sample
    xz ${long1%.xz} &
    xz ${long2%.xz}
    wait
else
    printf "Using existing $long1 and $long2.\n"
fi

short1=${sample}_1-250k.fastq.xz
short2=${sample}_2-250k.fastq.xz

if [ ! -e $short1 ]; then
    printf "Compressing...\n"
    xzcat -3 $long1 | head -1000000 | xz > $short1 &
    xzcat -3 $long2 | head -1000000 | xz > $short2
    wait
else
    printf "Using existing $short1 and $short2.\n"
fi

