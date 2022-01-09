#!/bin/sh -e

# Unfortunately, this sample looks like it's already been trimmed.
# Only 138 Nextera adapters found in almost a million reads and most not
# near the 3' end.  Probably natural sequences.  Find a dirty, publicly
# available sample for use here.
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
