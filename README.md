# Fastq-trim

## Description

Near-optimal read trimmer based on [biolibc](https://github.com/auerlab/biolibc)

The goal of fastq-trim is a read trimmer that runs in a fraction of the time
of popular trimmers and produces comparable results.
At present, fastq-trim uses simple alignment algorithms suitable for
typical analyses such as RNA-Seq or ATAC-Seq, where a small amount of
residual adapter content will not impact the downstream analysis.

The default is a simple function with two parameters, a minimum number of
bases matched and a maximum percentage of mismatches.
An exact-match function is also currently available, which runs slightly
faster and misses slightly more adapters.

These algorithms are not suitable for analyses that are highly sensitive to
adapter contamination.
However, fastq-trim is designed to allow dropping in a variety of alignment
functions.  More sophisticated functions can easily be added,
so fastq-trim can be easily adapted to perform trimming
for just about any purpose.

## Status

The current version has been tested on RNA-Seq and ATAC-Seq data.
The results so far are encouraging, with significantly better
performance than cutadapt or trimmomatic and results nearly identical to
cutadapt.  The default "smart" algorithm is the same one used by
cutadapt, except that it does not look for insertions and deletions (indels).

Diffing fastq-trim and cutadapt results reveals only
84 differences after trimming 1,000,000 reads.

Output from Test/compare-results.sh:

FastQ-Trim vs cutadapt:

```
Sequences in both files that differ:        68
Sequences only in fastq-trim output:        12
Sequences only in cutadapt output:           4
```

The differences are due mainly to cutadapt matching sequences with
indels (insertions
or deletions), while fastq-trim's default algorithm does not.  E.g.,
using an adapter sequence of CTGTCTCTTATA, one of the differences in
the trimmed reads is as follows:

```
fastq-trim: CCCCT ... CGCC CTGTCTCTTTATA CACATCTCC
cutadapt:   CCCCT ... CGCC
```

In this case, cutadapt assumed that CTGTCTCTTTATA was an adapter with
an inserted T, while fastq-trim assumed it is a natural sequence.
Whether this is more likely an adapter with an insertion or a natural
sequence similar to the adapter is up for debate.  With only 84 such
instances in 1,000,000 reads, it won't make any discernable difference
to a downstream RNA-Seq or ATAC-Seq analysis, where trimming isn't
even technically necessary (Liao, 2020, doi: 10.1093/nargab/lqaa068).

More functionality such as additional alignment algorithms and
command-line options will be added at a later date as time permits
and needs dictate.  Feel free to open an issue to request a new feature.

Resident memory (actual RAM) use peaks at around 2 MiB, which means
fastq-trim runs
almost entirely in cache RAM.  This is likely part of the reason
for its performance advantage over other tools.

Fastq-trim is currently single-threaded, as using additional cores will
not improve performance of the basic features.  Additional cores, if
available, may
be used by compression and decompression tools reading input files and
writing output files.

Run time with default parameters is only slightly longer than
```xzcat infile.xz | gzip > outfile.gz```, and fastq-trim on uncompressed
input and output actually outruns xzcat and gzip -1.  Multi-threading will
be considered if and when more computationally expensive features are added.

Results from running "./test.sh big" (in the Test directory) on a 2.9 GHz i5
are shown below.

```
Stats collected on an i5 2.9GHz 2-core, 4-hyperthread.

Peak CPU (including xzcat, gzip -1, pigz), wall time,
and peak memory (MiB, application only):

		CPU     Wall    Virtual Resident
Fastq-trim      260%    5.42    13      2
Cutadapt 2-core 310%    13.30   174     120
Trimmomatic     150%    21.03   3473    740
Cutadapt 1-core 140%    22.37   46      30

Detailed output:

Trimming with uncompressed input and output...

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Smart
  Maximum mismatch:  10%
  Filename:          temp-infile1.fastq
  Mode:              Single
  Adapter:           CTGTCTCTTATA

Read: 1000000  Adapter: 86768  Q < 20: 572761  Len < 30: 13012
	4.52 real         4.43 user         0.08 sys

All remaining tests use compressed input and output...

Timing compressed read and write without trimming...
	4.78 real         4.68 user         0.07 sys
	4.79 real         4.75 user         0.03 sys

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Exact
  Filename:          SRR1972918_1.fastq.xz
  Mode:              Single
  Adapter:           CTGTCTCTTATA

Read: 1000000  Adapter: 86180  Q < 20: 572761  Len < 30: 13000
	5.35 real        12.35 user         0.49 sys

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

Read: 1000000  Adapter: 86768  Q < 20: 572761  Len < 30: 13012
	5.42 real        13.88 user         0.37 sys

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

Read: 1000000  Adapter: 108722  Q < 20: 572761  Len < 30: 13346
	5.63 real        14.42 user         0.31 sys

Cutadapt 1 core...
status  in_reads        in_bp   too_short       too_long        too_many_n     out_reads        w/adapters      qualtrim_bp     out_bp
OK      1000000 101000000       13020   0       0       986980  86845   1005881388836129
       22.37 real        30.40 user         0.30 sys

Cutadapt 2 core...
status  in_reads        in_bp   too_short       too_long        too_many_n     out_reads        w/adapters      qualtrim_bp     out_bp
OK      1000000 101000000       13020   0       0       986980  86845   1005881388836129
       13.70 real        41.51 user         0.41 sys

TrimmomaticSE: Started with arguments:
 /dev/stdin SRR1972918_1-trimmed-trimmomatic.fastq.gz ILLUMINACLIP:nextera.fa:2:30:5 TRAILING:20 MINLEN:30
Automatically using 4 threads
Using Short Clipping Sequence: 'CTGTCTCTTATA'
ILLUMINACLIP: Using 0 prefix pairs, 1 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Reads: 1000000 Surviving: 987900 (98.79%) Dropped: 12100 (1.21%)
TrimmomaticSE: Completed successfully
       21.03 real        32.94 user         1.18 sys

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Smart
  Maximum mismatch:  10%
  Filename:          SRR1972918_2.fastq.xz
  Mode:              Single
  Adapter:           CTGTCTCTTATA

Read: 1000000  Adapter: 61965  Q < 20: 828265  Len < 30: 234054
	5.26 real        11.34 user         0.36 sys

Cutadapt 2 core reverse read...
status  in_reads        in_bp   too_short       too_long        too_many_n     out_reads        w/adapters      qualtrim_bp     out_bp
OK      1000000 101000000       234057  0       0       765943  62032   4028031456652196
       11.87 real        34.79 user         0.37 sys

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

Read: 2000000  Adapter: 97978  Q < 20: 1401026  Len < 30: 242967
	8.87 real        28.09 user         1.20 sys

Scanning for a portion of the adapter to flag any adapters missed due to
base substitutions.

Also scanning for random sequence TCGAACGGC as a control.

Raw data CTGTCTCTT   :                 63410
Raw data TCGAACGGC   :                  1525
Exact match output CTGTCTCTT   :         284
Exact match output TCGAACGGC   :        1304
Smart match 10 output CTGTCTCTT   :      185
Smart match 10 output TCGAACGGC   :     1304
Smart match 20 output CTGTCTCTT   :      103
Smart match 20 output TCGAACGGC   :     1300
Cutadapt output CTGTCTCTT   :            169
Cutadapt output TCGAACGGC   :           1304
Trimmomatic output CTGTCTCTT   :         103
Trimmomatic output TCGAACGGC   :        1308
```

## Design and Implementation

The code is organized following basic object-oriented design principals, but
implemented in C to minimize overhead and keep the source code accessible to
scientists who don't have time to master the complexities of C++.

Structures are treated as classes, with accessor macros and mutator functions
provided, so dependent applications and libraries need not access
structure members directly.  Since the C language cannot enforce this, it's
up to application programmers to exercise self-discipline.

For detailed coding standards, see
https://github.com/outpaddling/Coding-Standards/.

## Building and installing

Fastq-trim is intended to build cleanly in any POSIX environment on any CPU
architecture.  Please don't hesitate to open an issue if you encounter
problems on any Unix-like system.

Primary development is done on FreeBSD with clang, but the code is frequently
tested on Linux, MacOS, NetBSD, and OpenIndiana as well.  MS Windows is not supported,
unless using a POSIX environment such as Cygwin or Windows Subsystem for Linux.

The Makefile is designed to be friendly to package managers, such as
[Debian packages](https://www.debian.org/distrib/packages),
[FreeBSD ports](https://www.freebsd.org/ports/),
[MacPorts](https://www.macports.org/), [pkgsrc](https://pkgsrc.org/), etc.
End users should install using a package manager.  Note that pkgsrc can be used by anyone, on virtually any POSIX operating system, with or without administrator privileges..

I maintain a FreeBSD port and a pkgsrc package, which is sufficient to install
cleanly on virtually any POSIX platform.  If you would like to see a
package in another package manager, please consider creating a package
yourself.  This will be one of the easiest packages in the collection and
hence a good vehicle to learn how to create packages.

For an overview of available package managers, see the
[Repology website](https://repology.org/).

### Installing Fastq-trim on FreeBSD:

FreeBSD is a highly underrated platform for scientific computing, with over
2,000 scientific libraries and applications in the FreeBSD ports collection
(of more than 30,000 total), modern clang compiler, fully-integrated ZFS
filesystem, and renowned security, performance, and reliability.
FreeBSD has a somewhat well-earned reputation for being difficult to set up
and manage compared to user-friendly systems like [Ubuntu](https://ubuntu.com/).
However, if you're a little bit Unix-savvy, you can very quickly set up a
workstation, laptop, or VM using
[desktop-installer](http://www.acadix.biz/desktop-installer.php).
[GhostBSD](https://ghostbsd.org/) offers an experience very similar
to Ubuntu, but is build on FreeBSD rather than Debian Linux.  GhostBSD
packages lag behind FreeBSD by a few months, but this is not generally
an issue and there are workarounds.

To install the binary package on FreeBSD:

```
pkg install fastq-trim
```

You can just as easily build and install from source.  This is useful for
FreeBSD ports with special build options, for building with non-portable
optimizations such as -march=native, and for 
[work-in-progress ports](https://github.com/outpaddling/freebsd-ports-wip),
for which binary packages are not yet maintained.

```
cd /usr/ports/wip/fastq-trim && env CFLAGS='-march=native -O2' make install
```

### Installing via pkgsrc

pkgsrc is a cross-platform package manager that works on any Unix-like
platform. It is native to [NetBSD](https://www.netbsd.org/) and well-supported
on [Illumos](https://illumos.org/), [MacOS](https://www.apple.com/macos/),
[RHEL](https://www.redhat.com)/[CentOS](https://www.centos.org/), and
many other Linux distributions.
Using pkgsrc does not require admin privileges.  You can install a pkgsrc
tree in any directory to which you have write access and easily install any
of the nearly 20,000 packages in the collection.

The
[auto-pkgsrc-setup](https://github.com/outpaddling/auto-admin/blob/master/User-scripts/auto-pkgsrc-setup)
script will help you install pkgsrc in about 10 minutes.  Just download it
and run

```
sh auto-pkgsrc-setup
```

Then, assuming you selected current packages and the default prefix

```
source ~/Pkgsrc/pkg/etc/pkgsrc.sh   # Or pkgsrc.csh for csh or tcsh
cd ~/Pkgsrc/biology/fastq-trim
sbmake install clean clean-depends
```

See the pkgsrc documentation for more information.

Community support for pkgsrc is available through the
[pkgsrc-users](http://netbsd.org/mailinglists) mailing list.

### Building Fastq-trim locally

Below are caveman install instructions for development purposes, not
recommended for regular use.
Fastq-trim depends on [biolibc](https://github.com/auerlab/biolibc).
Install biolibc before attempting to build Fastq-trim.

1. Clone the repository
2. Run "make depend" to update Makefile.depend
3. Run "make install"

The default install prefix is ../local.  Clone Fastq-trim, biolibc and dependent
apps into sibling directories so that ../local represents a common path to all
of them.

To facilitate incorporation into package managers, the Makefile respects
standard make/environment variables such as CC, CFLAGS, PREFIX, etc.  

The library, headers, and man pages are installed under
`${DESTDIR}${PREFIX}`.  DESTDIR is empty by default and is primarily used by
package managers to stage installations.  PREFIX defaults to ../local.

Add-on libraries required for the build, such as biolibc, should be found
under either `${PREFIX}` or `${LOCALBASE}`, which defaults to `${PREFIX}`.
LOCALBASE can be set independently if you want to use libraries installed
by FreeBSD ports (/usr/local), MacPorts (/opt/local), pkgsrc (/usr/pkg), etc.

To install directly to /myprefix, assuming biolibc is installed there as well,
using a make variable:

```
make PREFIX=/myprefix clean depend install
```

Using an environment variable:

```
# C-shell and derivatives
setenv PREFIX /myprefix
make clean depend install

# Bourne shell and derivatives
PREFIX=/myprefix
export PREFIX
make clean depend install
```

View the Makefile for full details.
