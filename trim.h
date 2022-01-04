#define ILLUMINA_UNIVERSAL      "AGATCGGAAGAG"
#define ILLUMINA_SMALL_RNA_3P   "TGGAATTCTCGG"
#define ILLUMINA_SMALL_RNA_5P   "GATCGTCGGACT"
#define NEXTERA                 "CTGTCTCTTATA"
#define SOLID                   "CGCCTTGGCCGT"
#define DEFAULT_ADAPTER         ILLUMINA_UNIVERSAL

typedef struct
{
    char        *infile1;
    char        *outfile1;
    char        *infile2;
    char        *outfile2;
    FILE        *instream1;
    FILE        *outstream1;
    FILE        *instream2;
    FILE        *outstream2;
    char        *adapter;
    size_t      min_length;
    size_t      min_match;
    unsigned    min_qual;
    bool        verbose;
}   trim_t;

#include "trim-rvs.h"
#include "trim-accessors.h"
#include "trim-mutators.h"

