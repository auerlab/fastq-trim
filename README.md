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

The current version is the culmination of a few days work starting
from a blank slate, so gauge your expectations accordingly.  Mature tools
like cutadapt and trimmomatic will be more robust and feature-rich.

However, the results so far are encouraging, with better speed
than cutadapt and nearly identical results (diffing fastq-trim and cutadapt
results revealed only a few differences after trimming 250k reads).  Some
basic statistics from a larger sample with 1 million reads run on a
2.9 GHz i5 are below.  Note that fastq-trim is currently single-threaded.

```
Cores actually used (including xzcat, gzip, pigz):

Fastq-trim      3
Cutadapt 1-core 2
Cutadapt 2-core 4

Peak memory use (application only, not compression tools):

		Virtual Resident
Fastq-trim      13      2
Cutadapt 1-core 46      30
Cutadapt 2-core 58      40      ( each of 3 processes )
Trimmomatic     3473    740

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Exact
  Filename:          SRR1972918_1.fastq.xz
  Mode:              Single
  Adapter:           CTGTCTCTTATA

Read: 1000000  Adapter: 85032  Poly-A: 0  Q < 20: 572761  Len < 30: 92544
	5.91 real        13.06 user         0.32 sys

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Smart
  Maximum mismatch:  10%
  Filename:          SRR1972918_1.fastq.xz
  Mode:              Single
  Adapter:           CTGTCTCTTATA

Read: 1000000  Adapter: 85607  Poly-A: 0  Q < 20: 572761  Len < 30: 92553
	5.93 real        14.35 user         0.32 sys

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Smart
  Maximum mismatch:  20%
  Filename:          SRR1972918_1.fastq.xz
  Mode:              Single
  Adapter:           CTGTCTCTTATA

Read: 1000000  Adapter: 105995  Poly-A: 0  Q < 20: 572761  Len < 30: 92760
	6.35 real        14.90 user         0.27 sys

Cutadapt 1 core...
status  in_reads        in_bp   too_short       too_long        too_many_n     out_reads        w/adapters      qualtrim_bp     out_bp
OK      1000000 101000000       13020   0       0       986980  86845   1005881388836129
       23.09 real        31.17 user         0.25 sys

Cutadapt 2 core...
status  in_reads        in_bp   too_short       too_long        too_many_n     out_reads        w/adapters      qualtrim_bp     out_bp
OK      1000000 101000000       13020   0       0       986980  86845   1005881388836129
       13.83 real        43.64 user         0.52 sys

TrimmomaticSE: Started with arguments:
 /dev/stdin SRR1972918_1-trimmed-trimmomatic.fastq.gz ILLUMINACLIP:nextera.fa:2:30:5 TRAILING:20 MINLEN:30
Automatically using 4 threads
Using Short Clipping Sequence: 'CTGTCTCTTATA'
ILLUMINACLIP: Using 0 prefix pairs, 1 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Reads: 1000000 Surviving: 987900 (98.79%) Dropped: 12100 (1.21%)
TrimmomaticSE: Completed successfully
       20.94 real        33.22 user         0.95 sys

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Smart
  Maximum mismatch:  10%
  Filename:          SRR1972918_1.fastq.xz
  Filename:          SRR1972918_2.fastq.xz
  Mode:              Paired
  Adapters:          CTGTCTCTTATA AGATCGGAAGAG

Read: 1000000  Adapter: 92062  Poly-A: 0  Q < 20: 1401026  Len < 30: 545696
	8.08 real        26.41 user         0.55 sys

Scanning for a portion of the adapter to flag any adapters missed due to
base substitutions.

Also scanning for random sequence TCGAACGGC as a control.

Raw data CTGTCTCTT   :                 63410
Raw data TCGAACGGC   :                  1525
Exact match output CTGTCTCTT   :         275
Exact match output TCGAACGGC   :        1229
Smart match 10 output CTGTCTCTT   :      177
Smart match 10 output TCGAACGGC   :     1229
Smart match 20 output CTGTCTCTT   :       95
Smart match 20 output TCGAACGGC   :     1225
Cutadapt output CTGTCTCTT   :            169
Cutadapt output TCGAACGGC   :           1304
Trimmomatic output CTGTCTCTT   :         103
Trimmomatic output TCGAACGGC   :        1308
```
