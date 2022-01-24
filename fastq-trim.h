
#define ILLUMINA_UNIVERSAL      "AGATCGGAAGAG"
#define ILLUMINA_SMALL_RNA_3P   "TGGAATTCTCGG"
#define ILLUMINA_SMALL_RNA_5P   "GATCGTCGGACT"
#define NEXTERA                 "CTGTCTCTTATA"
#define SOLID                   "CGCCTTGGCCGT"
#define DEFAULT_ADAPTER         ILLUMINA_UNIVERSAL

#ifndef _BIOLIBC_FASTQ_H_
#include <biolibc/fastq.h>
#endif

#ifndef _BIOLIBC_ALIGN_H_
#include <biolibc/align.h>
#endif

#ifndef __bool_true_false_are_defined
#include <stdbool.h>
#endif

// Must typedef this to avoid confusing auto-gen-get-set in the structure
typedef size_t (* fastq_trim_afp_t)(const bl_align_t *params,
    const char *big, size_t big_len, const char *little, size_t little_len);

typedef struct
{
    bool        verbose;
    fastq_trim_afp_t  adapter_match_function;
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
    size_t      min_match;              // aggs:1:20
    size_t      polya_min_len;          // aggs:0:20
    unsigned    max_mismatch_percent;   // aggs:1:100
    unsigned    min_qual;               // aggs:1:60
    unsigned    phred_base;             // aggs:32:127
}   fastq_trim_t;

/* trim.c */
int     fastq_trim_single_reads(fastq_trim_t *tp);
int     fastq_trim_paired_reads(fastq_trim_t *tp);
void    fastq_trim_init(fastq_trim_t *tp);
int     fastq_trim_open_files(fastq_trim_t *tp, int arg, int argc, char *argv[]);
void    fastq_trim_close_files(fastq_trim_t *tp);

#include "fastq-trim-rvs.h"
#include "fastq-trim-accessors.h"
#include "fastq-trim-mutators.h"

