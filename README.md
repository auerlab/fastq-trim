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
2.9 GHz i5:

```
*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Exact
  Mode:              Single
  Adapter:           AGATCGGAAGAGCACAC

Reads: 32455547  Adapters: 990240  Qual < 20: 1604738  Len < 30: 620
      120.63 real       320.68 user         8.35 sys

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Smart
  Maximum mismatch:  10%
  Mode:              Single
  Adapter:           AGATCGGAAGAGCACAC

Reads: 32455547  Adapters: 993440  Qual < 20: 1604738  Len < 30: 629
      132.12 real       353.94 user        10.38 sys

Scanning for a portion of the adapter to flag any adapters missed due to
base substitutions.

Also scanning for random sequence of the same length for comparison.
Raw data AGATCGGAAGAG:                 39832
Raw data random:                         280
Exact match output AGATCGGAAGAG:         193
Exact match output random:               280
Smart match output AGATCGGAAGAG:          75
Smart match output random:               280
```
