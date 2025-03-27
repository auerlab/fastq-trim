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
#include <biolibc/align.h>
#include "fastq-trim.h"

// Explicit inlining makes no difference
size_t bl_fastq_find_polya_tail(bl_fastq_t *rec)

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

int     fastq_trim_single_reads(fastq_trim_t *tp)

{
    unsigned long   read_count,
		    adapter_count,
		    polya_count,
		    short_count,
		    low_qual_read_count,
		    low_qual_base_count,
		    adapter_pos_sum,
		    base_count;
    size_t          index, adapter_len;
    bl_fastq_t      rec;
    bl_align_t      align_params;
    
    fputs("  Read-type:         Single\n", stdout);
    fprintf(stdout,
	  "  Adapter:           %s\n\n", FASTQ_TRIM_ADAPTER1(tp));
    bl_fastq_init(&rec);

    adapter_len = strlen(tp->adapter1);
    bl_align_set_min_match(&align_params, tp->min_match);
    bl_align_set_max_mismatch_percent(&align_params, tp->max_mismatch_percent);
    
    read_count = adapter_count = polya_count = short_count = 
    low_qual_read_count = low_qual_base_count =
    adapter_pos_sum = base_count = 0;
    while ( bl_fastq_read(&rec, tp->instream1) == BL_READ_OK )
    {
	// Before trimming
	base_count += BL_FASTQ_SEQ_LEN(&rec);
	
	// Trim low-quality bases before adapters
	index = bl_fastq_find_3p_low_qual(&rec, tp->min_qual, tp->phred_base);
	if ( index != BL_FASTQ_QUAL_LEN(&rec) )
	{
	    low_qual_base_count += BL_FASTQ_QUAL_LEN(&rec) - index;
	    ++low_qual_read_count;
	    if ( tp->verbose )
		fprintf(stdout, "Low qual %s %s\n",
		    BL_FASTQ_SEQ(&rec) + index, BL_FASTQ_QUAL(&rec) + index);
	    bl_fastq_3p_trim(&rec, index);
	}

	index = bl_fastq_find_5p_low_qual(&rec, tp->min_qual, tp->phred_base);
	if ( index != -1 )
	{
	    low_qual_base_count += index + 1;
	    ++low_qual_read_count;
	    if ( tp->verbose )
	    {
		char    seq[index + 2], qual[index + 2];
		memcpy(seq, BL_FASTQ_SEQ(&rec), index + 1);
		memcpy(qual, BL_FASTQ_QUAL(&rec), index + 1);
		seq[index + 1] = qual[index + 1] = '\0';
		fprintf(stdout, "Low qual 5' %s %s\n", seq, qual);
	    }
	    bl_fastq_5p_trim(&rec, index);
	}

	index = tp->adapter_match_function(&align_params,
		    BL_FASTQ_SEQ(&rec), BL_FASTQ_SEQ_LEN(&rec),
		    tp->adapter1, adapter_len);
	if ( index != BL_FASTQ_SEQ_LEN(&rec) )
	{
	    adapter_pos_sum += index;
	    ++adapter_count;
	    if ( tp->verbose )
		fprintf(stdout, "Adapter  %s\n", BL_FASTQ_SEQ(&rec) + index);
	    bl_fastq_3p_trim(&rec, index);
	}

	if ( tp->polya_min_len != 0 )
	{
	    index = bl_fastq_find_polya_tail(&rec);
	    // Using unsigned and len could be < polya_min, so don't subtract
	    if ( index + tp->polya_min_len < BL_FASTQ_SEQ_LEN(&rec) )
	    {
		++polya_count;
		if ( tp->verbose )
		    fprintf(stdout, "Poly-A   %s\n",
			BL_FASTQ_SEQ(&rec) + index);
		bl_fastq_3p_trim(&rec, index);
	    }
	}

	if ( BL_FASTQ_SEQ_LEN(&rec) >= tp->min_length )
	    bl_fastq_write(&rec, tp->outstream1, BL_FASTQ_LINE_UNLIMITED);
	else
	{
	    if ( tp->verbose )
		fprintf(stdout, "Short    %zu %s\n",
			BL_FASTQ_SEQ_LEN(&rec), BL_FASTQ_SEQ(&rec));
	    ++short_count;
	}
	
	++read_count;
	
	// Display periodic progress only if running interactively
	if ( ! tp->verbose && (read_count % 100000 == 0) &&
	     isatty(fileno(stdout)) )
	{
	    fprintf(stdout,
		    "Read: %lu  Adapter: %lu  Poly-A: %lu  Q < %u: %lu  Len < %zu: %lu\r",
		    read_count, adapter_count, polya_count,
		    tp->min_qual, low_qual_read_count, tp->min_length, short_count);
	}
    }

    // Final results
    if ( read_count == 0 )
	printf("Input stream was empty.\n");
    else
	printf(
	    "\n\nReads:                             %10lu\n"
	    "Reads with adapters:               %10lu (%lu%%)\n"
	    "Reads with Poly-As:                %10lu (%lu%%)\n"
	    "Bases with Q < %u:                 %10lu (%lu%%)\n"
	    "Reads with low Q bases removed:    %10lu (%lu%%)\n"
	    "Reads < %zu bases after trimming:   %10lu (%lu%%)\n"
	    "Mean adapter position:             %10lu\n"
	    "Mean read length:                  %10lu\n",
	    read_count,
	    adapter_count, (adapter_count * 100 / read_count),
	    polya_count, (polya_count * 100 / read_count),
	    tp->min_qual, low_qual_base_count, (low_qual_base_count * 100 / base_count),
	    low_qual_read_count, (low_qual_read_count * 100 / read_count),
	    tp->min_length, short_count, (short_count * 100) / read_count,
	    adapter_pos_sum / adapter_count, base_count / read_count);
    bl_fastq_free(&rec);
    return EX_OK;
}


/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2022-01-03  Jason Bacon Begin
 ***************************************************************************/

int     fastq_trim_paired_reads(fastq_trim_t *tp)

{
    unsigned long   read_count,
		    adapter_count,
		    polya_count,
		    short_count,
		    low_qual_read_count,
		    low_qual_base_count,
		    adapter_pos_sum,
		    base_count;
    size_t          index, adapter_len[2];
    bl_fastq_t      rec[2];
    bl_align_t      align_params;
    int             s1, s2, c;
    char            *adapter[2];
    
    fputs("  Read-type:         Paired\n", stdout);
    fprintf(stdout,
	  "  Adapters:          %s %s\n\n", tp->adapter1, tp->adapter2);
    bl_fastq_init(&rec[0]);
    bl_fastq_init(&rec[1]);

    // Select adapters with loop index
    adapter[0] = tp->adapter1;
    adapter[1] = tp->adapter2;
    adapter_len[0] = strlen(adapter[0]);
    adapter_len[1] = strlen(adapter[1]);
    
    read_count = adapter_count = polya_count = short_count =
    low_qual_read_count = low_qual_base_count =
    adapter_pos_sum = base_count = 0;

    bl_align_set_min_match(&align_params, tp->min_match);
    bl_align_set_max_mismatch_percent(&align_params, tp->max_mismatch_percent);

    // Read from both files every iteration and break on error
    while ( true )
    {
	// FIXME: Explore using 2 threads here
	
	s1 = bl_fastq_read(&rec[0], tp->instream1);
	s2 = bl_fastq_read(&rec[1], tp->instream2);
	if ( (s1 != BL_READ_OK) || (s2 != BL_READ_OK) )
	    break;
	
	// Compare read names just for sanity checking
	if ( bl_fastq_name_cmp(&rec[0], &rec[1]) != 0 )
	{
	    fprintf(stderr, "fastq-trim: Paired files out of sync.\n");
	    fastq_trim_close_files(tp);
	    exit(EX_DATAERR);
	}
	
	for (c = 0; c <= 1; ++c)
	{
	    // Before trimming
	    base_count += BL_FASTQ_SEQ_LEN(&rec[c]);
	    
	    // Trim low quality bases before adapters
	    index = bl_fastq_find_3p_low_qual(&rec[c], tp->min_qual,
					  tp->phred_base);
	    if ( index != BL_FASTQ_QUAL_LEN(&rec[c]) )
	    {
		low_qual_base_count += BL_FASTQ_QUAL_LEN(&rec[c]) - index;
		++low_qual_read_count;
		if ( tp->verbose )
		    fprintf(stdout, "Low qual 3' %s %s\n",
			BL_FASTQ_SEQ(&rec[c]) + index,
			BL_FASTQ_QUAL(&rec[c]) + index);
		bl_fastq_3p_trim(&rec[c], index);
	    }

	    index = bl_fastq_find_5p_low_qual(&rec[c], tp->min_qual,
					  tp->phred_base);
	    if ( index != -1 )
	    {
		low_qual_base_count += index + 1;
		++low_qual_read_count;
		if ( tp->verbose )
		{
		    char    seq[index + 2], qual[index + 2];
		    memcpy(seq, BL_FASTQ_SEQ(&rec[c]), index + 1);
		    memcpy(qual, BL_FASTQ_QUAL(&rec[c]), index + 1);
		    seq[index + 1] = qual[index + 1] = '\0';
		    fprintf(stdout, "Low qual 5' %s %s\n", seq, qual);
		}
		bl_fastq_5p_trim(&rec[c], index);
	    }
    
	    index = tp->adapter_match_function(&align_params,
			BL_FASTQ_SEQ(&rec[c]), BL_FASTQ_SEQ_LEN(&rec[c]),
			adapter[c], adapter_len[c]);
	    if ( index != BL_FASTQ_SEQ_LEN(&rec[c]) )
	    {
		++adapter_count;
		adapter_pos_sum += index;
		if ( tp->verbose )
		    fprintf(stdout, "Adapter  %s\n",
			BL_FASTQ_SEQ(&rec[c]) + index);
		bl_fastq_3p_trim(&rec[c], index);
	    }

	    if ( tp->polya_min_len != 0 )
	    {
		index = bl_fastq_find_polya_tail(&rec[c]);
		// Using unsigned and len could be < polya_min, so don't subtract
		if ( index + tp->polya_min_len < BL_FASTQ_SEQ_LEN(&rec[c]) )
		{
		    ++polya_count;
		    if ( tp->verbose )
			fprintf(stdout, "Poly-A   %s\n",
			    BL_FASTQ_SEQ(&rec[c]) + index);
		    bl_fastq_3p_trim(&rec[c], index);
		}
	    }
	}        
	read_count += 2;

	/*
	 *  If either read is short, drop the pair.  Paired reads must be
	 *  kept in sync across the R1 and R2 files.
	 */
	if ( (BL_FASTQ_SEQ_LEN(&rec[0]) >= tp->min_length) &&
	     (BL_FASTQ_SEQ_LEN(&rec[1]) >= tp->min_length) )
	{
	    bl_fastq_write(&rec[0], tp->outstream1, BL_FASTQ_LINE_UNLIMITED);
	    bl_fastq_write(&rec[1], tp->outstream2, BL_FASTQ_LINE_UNLIMITED);
	}
	else
	{
	    if ( tp->verbose )
		fprintf(stdout, "Short    %zu %s\n"
				"         %zu %s\n",
			BL_FASTQ_SEQ_LEN(&rec[0]),
			BL_FASTQ_SEQ(&rec[1]),
			BL_FASTQ_SEQ_LEN(&rec[0]),
			BL_FASTQ_SEQ(&rec[1]));
	    ++short_count;
	}
	
	// Display periodic progress only if running interactively
	if ( ! tp->verbose && (read_count % 100000 == 0) &&
	     isatty(fileno(stdout)) )
	{
	    fprintf(stdout,
		    "Read: %lu  Adapter: %lu  Poly-A: %lu  Q < %u: %lu  Len < %zu: %lu\r",
		    read_count, adapter_count, polya_count,
		    tp->min_qual, low_qual_read_count, tp->min_length, short_count);
	}
    }
    
    // Final results
    fprintf(stdout,
	"\n\nReads:                             %10lu\n"
	"Reads with adapters:               %10lu (%lu%%)\n"
	"Reads with Poly-As:                %10lu (%lu%%)\n"
	"Bases with Q < %u:                 %10lu (%lu%%)\n"
	"Reads with low Q bases removed:    %10lu (%lu%%)\n"
	"Reads < %zu bases after trimming:   %10lu (%lu%%)\n"
	"Mean adapter position:             %10lu\n"
	"Mean read length:                  %10lu\n",
	read_count,
	adapter_count, (adapter_count * 100 / read_count),
	polya_count, (polya_count * 100 / read_count),
	tp->min_qual, low_qual_base_count, (low_qual_base_count * 100 / base_count),
	low_qual_read_count, (low_qual_read_count * 100 / read_count),
	tp->min_length, short_count, (short_count * 100) / read_count,
	adapter_pos_sum / adapter_count, base_count / read_count);

    bl_fastq_free(&rec[0]);
    bl_fastq_free(&rec[1]);
    return EX_OK;
}


/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2022-01-03  Jason Bacon Begin
 ***************************************************************************/

int     fastq_trim_open_files(fastq_trim_t *tp, int arg, int argc, char *argv[])

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
	    fastq_trim_close_files(tp);
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
	    fastq_trim_close_files(tp);
	    return EX_NOINPUT;
	}
	if ( arg == argc - 1 )
	{
	    tp->outfile2 = argv[arg++];
	    if ( (tp->outstream2 = xt_fopen(tp->outfile2, "w")) == NULL )
	    {
		fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0],
			tp->infile2, strerror(errno));
		fastq_trim_close_files(tp);
		return EX_CANTCREAT;
	    }
	    // paired_reads();
	}
	else
	    return EX_USAGE;
    }
    return EX_OK;
}


void    fastq_trim_close_files(fastq_trim_t *tp)

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


void    fastq_trim_init(fastq_trim_t *tp)

{
    tp->verbose = false;
    tp->adapter_match_function = bl_align_map_seq_sub;
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

