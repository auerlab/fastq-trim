# Fastq-trim
A lightening fast read trimmer.

This is a back-burner project to explore the possibility of developing a
significantly faster short read trimmer.
The ultimate goal is a trimmer that runs in a fraction of the time of
popular trimmers and produces good enough results so that a more
sophisticated trimmer won't make a meaningful difference to the downstream
analysis.

There is no sense of urgency since there are multiple highly-evolved trimmers
available that are fast enough for most purposes.  However, trimming can
take long enough to discourage experimenting, especially for those
who don't have access to an HPC cluster.

The current version is the culmination of roughly 2 days work starting
from a blank slate, so gauge your expectations accordingly.  Mature tools
like cutadapt and trimmomatic are much more robust and feature-rich.

However, the results so far are encouraging, with better speed
than cutadapt and nearly identical results (diffing fastq-trim and cutadapt
results revealed only a few differences after trimming 250k reads).  Some
basic statistics from a larger sample with ~32 million reads run on a
2.9 GHz i5 are below.  Note that fastq-trim is currently single-threaded.

```
Cores actually utilized (including xzcat, gzip, pigz):

Fastq-trim      3
Cutadapt 1-core 2
Cutadapt 2-core 4

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Exact
  Filename:          2M-R1.fastq.xz
  Mode:              Single
  Adapter:           AGATCGGAAGAGC

Reads: 2000000  Adapters: 101071  Qual < 20: 93395  Len < 30: 77
	6.35 real        16.99 user         0.91 sys

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Smart
  Maximum mismatch:  10%
  Filename:          2M-R1.fastq.xz
  Mode:              Single
  Adapter:           AGATCGGAAGAGC

Reads: 2000000  Adapters: 103096  Qual < 20: 93395  Len < 30: 89
	8.99 real        20.40 user         1.09 sys

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Smart
  Maximum mismatch:  20%
  Filename:          2M-R1.fastq.xz
  Mode:              Single
  Adapter:           AGATCGGAAGAGC

Reads: 2000000  Adapters: 142427  Qual < 20: 93395  Len < 30: 578
       10.10 real        21.34 user         1.30 sys

Cutadapt 1 core...
status  in_reads        in_bp   too_short       too_long        too_many_n      out_reads       w/adapters      qualtrim_bp     out_bp
OK      2000000 202000000       137     0       0       1999863 103341  180190  201023217
       44.10 real        60.33 user         0.72 sys

Cutadapt 2 core...
status  in_reads        in_bp   too_short       too_long        too_many_n      out_reads       w/adapters      qualtrim_bp     out_bp
OK      2000000 202000000       137     0       0       1999863 103341  180190  201023217
       22.49 real        61.15 user         1.48 sys

Scanning for a portion of the adapter to flag any adapters missed due to
base substitutions.

Also scanning for random sequence of the same length for comparison.

Raw data AGATCGGAAG  :                 23304
Raw data random:                         209
Exact match output AGATCGGAAG  :         208
Exact match output random:               209
Smart match 10 output AGATCGGAAG  :       83
Smart match 10 output random:            209
Smart match 20 output AGATCGGAAG  :       51
Smart match 20 output random:            208
Cutadapt output AGATCGGAAG  :             62
Cutadapt output random:                  209
```
