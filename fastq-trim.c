/***************************************************************************
 *  Description:
 *      Trim adapters and low-quality bases from FASTQ input
 *
 *  Arguments:
 *
 *  Returns:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include <xtend/file.h>
#include <biolibc/fastq.h>
#include "trim.h"

void    usage(char *argv[]);
size_t  bl_fastq_find_3p_adapter(const bl_fastq_t *read, const char *adapter, size_t min_match);
size_t  bl_fastq_find_3p_qual(const bl_fastq_t *read, unsigned min_qual, unsigned phred_base);
size_t  bl_fastq_3p_trim(bl_fastq_t *read, size_t new_len);
int     trim_single_reads(trim_t *tp);
int     trim_paired_reads(trim_t *tp);
void    trim_init(trim_t *tp);
int     trim_open_files(trim_t *tp, int arg, int argc, char *argv[]);
int     main(int argc,char *argv[])

{
    int     arg;
    char    *end;
    trim_t  tp;
    
    trim_init(&tp);
    
    if ( (argc == 2) && (strcmp(argv[1], "--help") == 0) )
	usage(argv);
    for (arg = 1; (arg < argc) && (*argv[arg] == '-'); ++arg)
    {
	if ( strcmp(argv[arg], "--verbose") == 0 )
	    trim_set_verbose(&tp, true);
	else if ( strcmp(argv[arg], "--3p-adapter") == 0 )
	{
	    // FIXME: tolower() adapter ahead of time and tolower() the read
	    // bases in bl_fastq_find_3p_adapter()
	    trim_set_adapter(&tp, argv[++arg]);
	}
	else if ( strcmp(argv[arg], "--min-match") == 0 )
	{
	    trim_set_min_match(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--min-qual") == 0 )
	{
	    trim_set_min_qual(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--min-length") == 0 )
	{
	    trim_set_min_length(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else
	    usage(argv);
    }
    
    if ( trim_open_files(&tp, arg, argc, argv) == EX_OK )
    {
	// FIXME: Support paired reads
	trim_single_reads(&tp);
    }
    
    // FIXME: trim_close_files()
    xt_fclose(tp.outstream1);
    xt_fclose(tp.instream1);
    return EX_OK;
}




/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2022-01-03  Jason Bacon Begin
 ***************************************************************************/

int     trim_open_files(trim_t *tp, int arg, int argc, char *argv[])

{
    // infile1 must be named before outfile1
    if ( arg < argc )
    {
	tp->infile1 = argv[arg++];
	if ( (tp->instream1 = xt_fopen(tp->infile1, "r")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0],
		    tp->infile1, strerror(errno));
	    return EX_NOINPUT;
	}
    }
    
    if ( arg < argc )
    {
	tp->outfile1 = argv[arg++];
	if ( (tp->outstream1 = xt_fopen(tp->outfile1, "w")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0],
		    tp->infile1, strerror(errno));
	    xt_fclose(tp->instream1);
	    return EX_CANTCREAT;
	}
    }
    
    // Infile2 and outfile2 must both be named or neither
    if ( arg < argc )
    {
	tp->infile2 = argv[arg++];
	if ( (tp->instream2 = xt_fopen(tp->infile2, "r")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0],
		    tp->infile2, strerror(errno));
	    return EX_NOINPUT;
	}
	if ( arg == argc - 1 )
	{
	    tp->outfile2 = argv[arg++];
	    if ( (tp->outstream2 = xt_fopen(tp->outfile2, "w")) == NULL )
	    {
		fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0],
			tp->infile2, strerror(errno));
		return EX_CANTCREAT;
	    }
	    // paired_reads();
	}
	else
	    usage(argv);
    }
    return EX_OK;
}


void    trim_init(trim_t *tp)

{
    tp->infile1 = NULL;
    tp->outfile1 = NULL;
    tp->infile2 = NULL;
    tp->outfile2 = NULL;
    tp->instream1 = stdin;
    tp->outstream1 = stdout;
    tp->instream2 = NULL;
    tp->outstream2 = NULL;
    tp->adapter = NULL;
    tp->min_length = 30;
    tp->min_match = 3;
    tp->min_qual = 20;
    tp->verbose = false;
}


/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2022-01-03  Jason Bacon Begin
 ***************************************************************************/

int     trim_single_reads(trim_t *tp)

{
    unsigned long   record_count,
		    adapter_count,
		    short_count,
		    low_qual_count;
    size_t          index;
    bl_fastq_t      fastq_rec1;
    
    bl_fastq_init(&fastq_rec1);
    record_count = adapter_count = short_count = low_qual_count = 0;
    while ( bl_fastq_read(tp->instream1, &fastq_rec1) == BL_READ_OK )
    {
	index = bl_fastq_find_3p_adapter(&fastq_rec1, tp->adapter, tp->min_match);
	if ( BL_FASTQ_SEQ_AE(&fastq_rec1, index) != '\0' )
	{
	    ++adapter_count;
	    if ( tp->verbose )
		fprintf(stderr, "Adapter  %s\n",
			BL_FASTQ_SEQ(&fastq_rec1) + index);
	    bl_fastq_3p_trim(&fastq_rec1, index);
	}
	
	// FIXME: Support other PHRED bases
	index = bl_fastq_find_3p_qual(&fastq_rec1, tp->min_qual, 33);
	if ( BL_FASTQ_SEQ_AE(&fastq_rec1, index) != '\0' )
	{
	    ++low_qual_count;
	    if ( tp->verbose )
		fprintf(stderr, "Low qual %s\n",
			BL_FASTQ_SEQ(&fastq_rec1) + index);
	    bl_fastq_3p_trim(&fastq_rec1, index);
	}
	
	//fprintf(stderr, "%zu\n", BL_FASTQ_SEQ_LEN(&fastq_rec1));
	if ( BL_FASTQ_SEQ_LEN(&fastq_rec1) >= tp->min_length )
	{
	    bl_fastq_write(tp->outstream1, &fastq_rec1, BL_FASTQ_SEQ_LEN(&fastq_rec1));
	}
	else
	{
	    if ( tp->verbose )
		fprintf(stderr, "Short    %zu %s\n",
			BL_FASTQ_SEQ_LEN(&fastq_rec1),
			BL_FASTQ_SEQ(&fastq_rec1));
	    ++short_count;
	}
	
	++record_count;
	if ( ! tp->verbose && (record_count % 100000 == 0) )
	{
	    fprintf(stderr,
		    "Reads: %lu  Adapters: %lu  Qual < %u: %lu  Len < %zu: %lu\r",
		    record_count, adapter_count, tp->min_qual,
		    low_qual_count, tp->min_length, short_count);
	}
    }
    fprintf(stderr,
	    "Reads: %lu  Adapters: %lu  Qual < %u: %lu  Len < %zu: %lu\n",
	    record_count, adapter_count, tp->min_qual, low_qual_count,
	    tp->min_length, short_count);
    bl_fastq_free(&fastq_rec1);
    return EX_OK;
}


/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2022-01-03  Jason Bacon Begin
 ***************************************************************************/

int     trim_paired_reads(trim_t *tp)

{
    unsigned long   record_count,
		    adapter_count,
		    short_count,
		    low_qual_count;
    size_t          index;
    bl_fastq_t      fastq_rec1;
    
    bl_fastq_init(&fastq_rec1);
    record_count = adapter_count = short_count = low_qual_count = 0;
    while ( bl_fastq_read(tp->instream1, &fastq_rec1) == BL_READ_OK )
    {
	index = bl_fastq_find_3p_adapter(&fastq_rec1, tp->adapter, tp->min_match);
	if ( BL_FASTQ_SEQ_AE(&fastq_rec1, index) != '\0' )
	{
	    ++adapter_count;
	    if ( tp->verbose )
		fprintf(stderr, "Adapter  %s\n",
			BL_FASTQ_SEQ(&fastq_rec1) + index);
	    bl_fastq_3p_trim(&fastq_rec1, index);
	}
	
	// FIXME: Support other PHRED bases
	index = bl_fastq_find_3p_qual(&fastq_rec1, tp->min_qual, 33);
	if ( BL_FASTQ_SEQ_AE(&fastq_rec1, index) != '\0' )
	{
	    ++low_qual_count;
	    if ( tp->verbose )
		fprintf(stderr, "Low qual %s\n",
			BL_FASTQ_SEQ(&fastq_rec1) + index);
	    bl_fastq_3p_trim(&fastq_rec1, index);
	}
	
	//fprintf(stderr, "%zu\n", BL_FASTQ_SEQ_LEN(&fastq_rec1));
	if ( BL_FASTQ_SEQ_LEN(&fastq_rec1) >= tp->min_length )
	{
	    bl_fastq_write(tp->outstream1, &fastq_rec1, BL_FASTQ_SEQ_LEN(&fastq_rec1));
	}
	else
	{
	    if ( tp->verbose )
		fprintf(stderr, "Short    %zu %s\n",
			BL_FASTQ_SEQ_LEN(&fastq_rec1),
			BL_FASTQ_SEQ(&fastq_rec1));
	    ++short_count;
	}
	
	++record_count;
	if ( ! tp->verbose && (record_count % 100000 == 0) )
	{
	    fprintf(stderr,
		    "Reads: %lu  Adapters: %lu  Qual < %u: %lu  Len < %zu: %lu\r",
		    record_count, adapter_count, tp->min_qual,
		    low_qual_count, tp->min_length, short_count);
	}
    }
    fprintf(stderr,
	    "Reads: %lu  Adapters: %lu  Qual < %u: %lu  Len < %zu: %lu\n",
	    record_count, adapter_count, tp->min_qual, low_qual_count,
	    tp->min_length, short_count);
    bl_fastq_free(&fastq_rec1);
    return EX_OK;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate start of a low-quality 3' end in a FASTQ read.
 *  
 *  Arguments:
 *      read        FASTQ read to be searched
 *      min_qual    Minimum quality of bases to keep
 *
 *  Returns:
 *      Index of first low-quality base at the 3' end if found,
 *      index of NULL terminator otherwise
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_find_3p_qual(const bl_fastq_t *read, unsigned min_qual,
			unsigned phred_base)

{
    const char  *p;
    size_t      index = read->qual_len;
    
    for (p = read->qual + read->qual_len - 1;
	 (p >= read->qual) && (toupper(*p) - phred_base <= min_qual); --p)
    {
	index = p - read->qual;
	// fprintf(stderr, "%c %u\n", read->seq[index], read->qual[index] - phred_base);
    }
    return index;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate adapter (or a portion thereof if at the end) starting at the
 *      3' end of the read.
 *  
 *  Arguments:
 *      read        FASTQ read to be searched
 *      adapter     Adapter sequence to be located
 *      min_match   Minimum number of characters to match in adapter
 *
 *  Returns:
 *      Index of adapter sequence if found, index of NULL terminator otherwise
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_find_3p_adapter(const bl_fastq_t *read, const char *adapter,
			 size_t min_match)

{
    const char *start, *pr, *pa;
    
    for (start = read->seq + read->seq_len - 1; start >= read->seq; --start)
    {
	for (pr = start, pa = adapter; toupper(*pr) == *pa; ++pr, ++pa)
	    ;
	if ( ((*pr == '\0') && (pr - start >= min_match)) || (*pa == '\0') )
	    return start - read->seq;
    }
    return read->seq_len;   // Location of '\0' terminator
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Trim the 3' end of a FASTQ read.
 *  
 *  Arguments:
 *      read        FASTQ read to be searched
 *      new_len     New length and location of the null terminators
 *
 *  Returns:
 *      BL_DATA_OK if new_len is between 0 and original length,
 *      BL_DATA_INVALID otherwise.
 *
 *  Examples:
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_3p_trim(bl_fastq_t *read, size_t new_len)

{
    if ( (new_len >= 0) && (new_len <= read->seq_len) )
    {
	read->seq_len = read->qual_len = new_len;
	read->seq[new_len] = read->qual[new_len] = '\0';
	// FIXME: realloc?
	return BL_DATA_OK;
    }
    else
	return BL_DATA_INVALID;
}


void    usage(char *argv[])

{
    fprintf(stderr,
	    "Usage:\n\n"
	    "%s --help\n"
	    "%s\n"
	    "    [--3p-adapter seq] [--min-qual N] [--min-length N]\n"
	    "    [infile1.fastq[.xz|.bz2|.gz]] [outfile1.fastq[.xz|.bz2|.gz]]\n\n"
	    "    [infile2.fastq[.xz|.bz2|.gz]] [outfile2.fastq[.xz|.bz2|.gz]]\n\n",
	    argv[0], argv[0]);
    exit(EX_USAGE);
}

