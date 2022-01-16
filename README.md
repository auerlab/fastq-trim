# Fastq-trim

## Description

Near-optimal read trimmer based on [biolibc](https://github.com/auerlab/biolibc)

This is a back-burner project to explore the possibility of developing a
significantly faster read trimmer.
The ultimate goal is a trimmer that runs in a fraction of the time of
popular trimmers and produces comparable results.

There is no sense of urgency since there are multiple highly-evolved trimmers
available that are fast enough for most purposes.  However, trimming can
take long enough to discourage experimenting with trim parameters,
especially for those who don't have access to an HPC cluster.

The current version is the culmination of a few days work starting
from a blank slate, so gauge your expectations accordingly.  Mature tools
like cutadapt and trimmomatic will be more robust and feature-rich.

However, the results so far are encouraging, with significantly better
speed than cutadapt or trimmomatic and nearly identical results (diffing
fastq-trim and cutadapt results revealed only a few differences after
trimming 250k reads).

Resident memory peaks at around 2 MiB, which means fastq-trim will run
almost entirely in cache RAM.

Run time with default parameters is only slightly longer than
```xzcat infile.xz | gzip > outfile.gz```, so there isn't much room for
improvement.

Results from running "./test.sh big" (in the Test directory) on a 2.9 GHz i5
are shown below.  Note that fastq-trim is currently single-threaded.

```
Stats collected on an i5 2.9GHz 2-core, 4-hyperthread.

Peak CPU (including xzcat, gzip -1, pigz), wall time, and peak memory (MiB):

		CPU     Wall    Virtual Resident
Fastq-trim      260%    5.42    13      2
Cutadapt 1-core 140%    22.37   46      30
Cutadapt 2-core 310%    13.30   174     120
Trimmomatic     150%    21.03   3473    740

Detailed output:

Timing read and write without trimming...
	4.80 real         4.77 user         0.00 sys
	4.82 real         4.75 user         0.06 sys

*** FASTQ TRIM ***

  Minimum match:     3
  Minimum quality:   20
  Minimum length:    30
  Phred base:        33
  Adapter matching:  Exact
  Filename:          SRR1972918_1.fastq.xz
  Mode:              Single
  Adapter:           CTGTCTCTTATA

Read: 1000000  Adapter: 86180  Poly-A: 0  Q < 20: 572761  Len < 30: 13000
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

Read: 1000000  Adapter: 86768  Poly-A: 0  Q < 20: 572761  Len < 30: 13012
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

Read: 1000000  Adapter: 108722  Poly-A: 0  Q < 20: 572761  Len < 30: 13346
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

Read: 1000000  Adapter: 61965  Poly-A: 0  Q < 20: 828265  Len < 30: 234054
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

Read: 2000000  Adapter: 97978  Poly-A: 0  Q < 20: 1401026  Len < 30: 242967
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
tested on CentOS, MacOS, and NetBSD as well.  MS Windows is not supported,
unless using a POSIX environment such as Cygwin or Windows Subsystem for Linux.

The Makefile is designed to be friendly to package managers, such as
[Debian packages](https://www.debian.org/distrib/packages),
[FreeBSD ports](https://www.freebsd.org/ports/),
[MacPorts](https://www.macports.org/), [pkgsrc](https://pkgsrc.org/), etc.
End users should install via one of these if at all possible.

I maintain a FreeBSD port and a pkgsrc package, which is sufficient to install
cleanly on virtually any POSIX platform.  If you would like to see a
package in another package manager, please consider creating a package
yourself.  This will be one of the easiest packages in the collection and
hence a good vehicle to learn how to create packages.

### Installing Fastq-trim on FreeBSD:

FreeBSD is a highly underrated platform for scientific computing, with over
1,900 scientific libraries and applications in the FreeBSD ports collection
(of more than 30,000 total), modern clang compiler, fully-integrated ZFS
filesystem, and renowned security, performance, and reliability.
FreeBSD has a somewhat well-earned reputation for being difficult to set up
and manage compared to user-friendly systems like [Ubuntu](https://ubuntu.com/).
However, if you're a little bit Unix-savvy, you can very quickly set up a
workstation, laptop, or VM using
[desktop-installer](http://www.acadix.biz/desktop-installer.php).

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
of the nearly 20,000 packages in the collection.  The
[auto-pkgsrc-setup](http://netbsd.org/~bacon/) script can assist you with
basic setup.

First bootstrap pkgsrc using auto-pkgsrc-setup or any
other method.  Then run the following commands:

```
cd pkgsrc-dir/wip/fastq-trim
bmake install clean
```

There may also be binary packages available for your platform.  If this is
the case, you can install by running:

```
pkgin install fastq-trim
```

See the [Joyent Cloud Services Site](https://pkgsrc.joyent.com/) for
available package sets.

### Building Fastq-trim locally

Below are cave man install instructions for development purposes, not
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

Add-on libraries required for the build, such as biolibc, should be found
under ${LOCALBASE}, which defaults to ../local.
The library, headers, and man pages are installed under
${DESTDIR}${PREFIX}.  DESTDIR is empty by default and is primarily used by
package managers to stage installations.  PREFIX defaults to ${LOCALBASE}.

To install directly to /myprefix, assuming biolibc is installed there as well,
using a make variable:

```
make LOCALBASE=/myprefix clean depend install
```

Using an environment variable:

```
# C-shell and derivatives
setenv LOCALBASE /myprefix
make clean depend install

# Bourne shell and derivatives
LOCALBASE=/myprefix
export LOCALBASE
make clean depend install
```

View the Makefile for full details.
