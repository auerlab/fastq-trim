
#define ILLUMINA_UNIVERSAL      "AGATCGGAAGAG"
#define ILLUMINA_SMALL_RNA_3P   "TGGAATTCTCGG"
#define ILLUMINA_SMALL_RNA_5P   "GATCGTCGGACT"
#define NEXTERA                 "CTGTCTCTTATA"
#define SOLID                   "CGCCTTGGCCGT"
#define DEFAULT_ADAPTER         ILLUMINA_UNIVERSAL

#ifndef _BIOLIBC_FASTQ_H_
#include <biolibc/fastq.h>
#endif

#ifndef __bool_true_false_are_defined
#include <stdbool.h>
#endif

// Must typedef this to avoid confusing auto-gen-get-set in the structure
typedef size_t (* trim_afp_t)(const bl_fastq_t *read,
	    const char *adapter, size_t min_match, unsigned max_mismatch_percent);

typedef struct
{
    bool        verbose;
    trim_afp_t  adapter_match_function;
    char        *infile1;
    char        *outfile1;
    char        *infile2;
    char        *outfile2;
    FILE        *instream1;
    FILE        *outstream1;
    FILE        *instream2;
    FILE        *outstream2;
    char        *adapter1;
    char        *adapter2;
    size_t      min_length;
    size_t      min_match;
    size_t      polya_min_len;
    unsigned    max_mismatch_percent;
    unsigned    min_qual;
    unsigned    phred_base;
}   trim_t;

/* trim.c */
int     trim_single_reads(trim_t *tp);
int     trim_paired_reads(trim_t *tp);
void    trim_init(trim_t *tp);
int     trim_open_files(trim_t *tp, int arg, int argc, char *argv[]);
void    trim_close_files(trim_t *tp);

#include "trim-rvs.h"
#include "trim-accessors.h"
#include "trim-mutators.h"

