#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>         // isatty()
#include <ctype.h>
#include <xtend/file.h>
#include <xtend/string.h>
#include <biolibc/fastq.h>
#include "trim.h"

static inline size_t bl_fastq_find_polya_tail(bl_fastq_t *rec)

{
    ssize_t c;
    
    for (c = rec->seq_len - 1; toupper(rec->seq[c] == 'A') && (c >= 0); --c)
	;
    return c + 1;
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
		    polya_count,
		    short_count,
		    low_qual_count;
    size_t          index;
    bl_fastq_t      fastq_rec;
    
    fputs("  Mode:              Single\n", stderr);
    fprintf(stderr,
	  "  Adapter:           %s\n\n", TRIM_ADAPTER1(tp));
    bl_fastq_init(&fastq_rec);
    record_count = adapter_count = polya_count = short_count = low_qual_count = 0;
    while ( bl_fastq_read(&fastq_rec, tp->instream1) == BL_READ_OK )
    {
	// Trim low-quality bases before adapters
	index = bl_fastq_find_3p_low_qual(&fastq_rec, tp->min_qual, tp->phred_base);
	if ( index != BL_FASTQ_QUAL_LEN(&fastq_rec) )
	{
	    ++low_qual_count;
	    if ( tp->verbose )
		fprintf(stderr, "Low qual %s %s\n",
		    BL_FASTQ_SEQ(&fastq_rec) + index,
		    BL_FASTQ_QUAL(&fastq_rec) + index);
	    bl_fastq_3p_trim(&fastq_rec, index);
	}

	index = tp->adapter_match_function(&fastq_rec, tp->adapter1,
		    tp->min_match, tp->max_mismatch_percent);
	if ( index != BL_FASTQ_SEQ_LEN(&fastq_rec) )
	{
	    ++adapter_count;
	    if ( tp->verbose )
		fprintf(stderr, "Adapter  %s\n",
		    BL_FASTQ_SEQ(&fastq_rec) + index);
	    bl_fastq_3p_trim(&fastq_rec, index);
	}

	if ( tp->polya_min_len != 0 )
	{
	    index = bl_fastq_find_polya_tail(&fastq_rec);
	    // Using unsigned and len could be < 10, so don't subtract
	    if ( index + tp->polya_min_len < BL_FASTQ_SEQ_LEN(&fastq_rec) )
	    {
		++polya_count;
		if ( tp->verbose )
		    fprintf(stderr, "Poly-A   %s\n",
			BL_FASTQ_SEQ(&fastq_rec) + index);
		bl_fastq_3p_trim(&fastq_rec, index);
	    }
	}

	if ( BL_FASTQ_SEQ_LEN(&fastq_rec) >= tp->min_length )
	    bl_fastq_write(&fastq_rec, tp->outstream1, BL_FASTQ_LINE_UNLIMITED);
	else
	{
	    if ( tp->verbose )
		fprintf(stderr, "Short    %zu %s\n",
			BL_FASTQ_SEQ_LEN(&fastq_rec),
			BL_FASTQ_SEQ(&fastq_rec));
	    ++short_count;
	}
	
	++record_count;
	
	// isatty() oddly returns true under SLURM
	if ( ! tp->verbose && (record_count % 100000 == 0) &&
	     isatty(fileno(stderr)) && (getenv("SLURM_JOB_ID") == NULL) )
	{
	    fprintf(stderr,
		    "Read: %lu  Adapter: %lu  Poly-A: %lu  Q < %u: %lu  Len < %zu: %lu\r",
		    record_count, adapter_count, polya_count,
		    tp->min_qual, low_qual_count, tp->min_length, short_count);
	}
    }
    fprintf(stderr,
	    "Read: %lu  Adapter: %lu  Poly-A: %lu  Q < %u: %lu  Len < %zu: %lu\n",
	    record_count, adapter_count, polya_count,
	    tp->min_qual, low_qual_count, tp->min_length, short_count);
    bl_fastq_free(&fastq_rec);
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
		    polya_count,
		    short_count,
		    low_qual_count;
    size_t          index;
    bl_fastq_t      fastq_rec[2];
    int             s1, s2, c;
    char            *adapter[2];
    
    fputs("  Mode:              Paired\n", stderr);
    fprintf(stderr,
	  "  Adapters:          %s %s\n\n", tp->adapter1, tp->adapter2);
    bl_fastq_init(&fastq_rec[0]);
    bl_fastq_init(&fastq_rec[1]);

    // Select adapters with loop index
    adapter[0] = tp->adapter1;
    adapter[1] = tp->adapter2;
    
    record_count = adapter_count = polya_count = short_count = low_qual_count = 0;

    // Read from both files every iteration and break on error
    while ( true )
    {
	// FIXME: Explore using 2 threads here
	
	s1 = bl_fastq_read(&fastq_rec[0], tp->instream1);
	s2 = bl_fastq_read(&fastq_rec[1], tp->instream2);
	if ( (s1 != BL_READ_OK) || (s2 != BL_READ_OK) )
	    break;
	
	// Compare read names just for sanity checking
	if ( bl_fastq_name_cmp(&fastq_rec[0], &fastq_rec[1]) != 0 )
	{
	    fprintf(stderr, "fastq-trim: Paired files out of sync.\n");
	    trim_close_files(tp);
	    exit(EX_DATAERR);
	}
	
	for (c = 0; c <= 1; ++c)
	{
	    // Trim low quality bases before adapters
	    index = bl_fastq_find_3p_low_qual(&fastq_rec[c], tp->min_qual,
					  tp->phred_base);
	    if ( index != BL_FASTQ_QUAL_LEN(&fastq_rec[c]) )
	    {
		++low_qual_count;
		if ( tp->verbose )
		    fprintf(stderr, "Low qual %s %s\n",
			BL_FASTQ_SEQ(&fastq_rec[c]) + index,
			BL_FASTQ_QUAL(&fastq_rec[c]) + index);
		bl_fastq_3p_trim(&fastq_rec[c], index);
	    }
    
	    index = tp->adapter_match_function(&fastq_rec[c], adapter[c],
			tp->min_match, tp->max_mismatch_percent);
	    if ( index != BL_FASTQ_SEQ_LEN(&fastq_rec[c]) )
	    {
		++adapter_count;
		if ( tp->verbose )
		    fprintf(stderr, "Adapter  %s\n",
			    BL_FASTQ_SEQ(&fastq_rec[c]) + index);
		bl_fastq_3p_trim(&fastq_rec[c], index);
	    }

	    if ( tp->polya_min_len != 0 )
	    {
		index = bl_fastq_find_polya_tail(&fastq_rec[c]);
		// Using unsigned and len could be < 10, so don't subtract
		if ( index + tp->polya_min_len < BL_FASTQ_SEQ_LEN(&fastq_rec[c]) )
		{
		    ++polya_count;
		    if ( tp->verbose )
			fprintf(stderr, "Poly-A   %s\n",
			    BL_FASTQ_SEQ(&fastq_rec[c]) + index);
		    bl_fastq_3p_trim(&fastq_rec[c], index);
		}
	    }
	}        

	/*
	 *  If either read is short, drop the pair.  Paired reads must be
	 *  kept in sync across the R1 and R2 files.
	 */
	if ( (BL_FASTQ_SEQ_LEN(&fastq_rec[0]) >= tp->min_length) &&
	     (BL_FASTQ_SEQ_LEN(&fastq_rec[1]) >= tp->min_length) )
	{
	    bl_fastq_write(&fastq_rec[0], tp->outstream1, BL_FASTQ_LINE_UNLIMITED);
	    bl_fastq_write(&fastq_rec[1], tp->outstream2, BL_FASTQ_LINE_UNLIMITED);
	}
	else
	{
	    if ( tp->verbose )
		fprintf(stderr, "Short    %zu %s\n"
				"         %zu %s\n",
			BL_FASTQ_SEQ_LEN(&fastq_rec[0]),
			BL_FASTQ_SEQ(&fastq_rec[1]),
			BL_FASTQ_SEQ_LEN(&fastq_rec[0]),
			BL_FASTQ_SEQ(&fastq_rec[1]));
	    ++short_count;
	}
	
	++record_count;
	if ( ! tp->verbose && (record_count % 100000 == 0) )
	{
	    fprintf(stderr,
		    "Read: %lu  Adapter: %lu  Poly-A: %lu  Q < %u: %lu  Len < %zu: %lu\r",
		    record_count, adapter_count, polya_count,
		    tp->min_qual, low_qual_count, tp->min_length, short_count);
	}
    }
    fprintf(stderr,
	    "Read: %lu  Adapter: %lu  Poly-A: %lu  Q < %u: %lu  Len < %zu: %lu\n",
	    record_count, adapter_count, polya_count,
	    tp->min_qual, low_qual_count, tp->min_length, short_count);

    bl_fastq_free(&fastq_rec[0]);
    bl_fastq_free(&fastq_rec[1]);
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
	    trim_close_files(tp);
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
	    trim_close_files(tp);
	    return EX_NOINPUT;
	}
	if ( arg == argc - 1 )
	{
	    tp->outfile2 = argv[arg++];
	    if ( (tp->outstream2 = xt_fopen(tp->outfile2, "w")) == NULL )
	    {
		fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0],
			tp->infile2, strerror(errno));
		trim_close_files(tp);
		return EX_CANTCREAT;
	    }
	    // paired_reads();
	}
	else
	    return EX_USAGE;
    }
    return EX_OK;
}


void    trim_close_files(trim_t *tp)

{
    if ( tp->infile1 != NULL )
	xt_fclose(tp->instream1);
    if ( tp->outfile1 != NULL )
	xt_fclose(tp->outstream1);
    if ( tp->infile2 != NULL )
	xt_fclose(tp->instream2);
    if ( tp->outfile2 != NULL )
	xt_fclose(tp->outstream2);
}


void    trim_init(trim_t *tp)

{
    tp->verbose = false;
    tp->adapter_match_function = bl_fastq_find_adapter_smart;
    tp->infile1 = NULL;
    tp->outfile1 = NULL;
    tp->infile2 = NULL;
    tp->outfile2 = NULL;
    tp->instream1 = stdin;
    tp->outstream1 = stdout;
    tp->instream2 = NULL;
    tp->outstream2 = NULL;
    tp->adapter1 = strdup(ILLUMINA_UNIVERSAL);
    tp->adapter2 = strdup(ILLUMINA_UNIVERSAL);
    tp->min_length = 30;
    tp->min_match = 3;
    tp->polya_min_len = 0;
    tp->max_mismatch_percent = 10;
    tp->min_qual = 20;
    tp->phred_base = 33;
}

