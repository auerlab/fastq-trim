/***************************************************************************
 *  Description:
 *      Trim adapters and low-quality bases from FASTQ input
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <xtend/string.h>
#include <biolibc/fastq.h>
#include "fastq-trim.h"

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    int     arg, status;
    char    *end;
    fastq_trim_t  tp;
    
    fastq_trim_init(&tp);
    
    if ( (argc < 2) || ((argc == 2) && (strcmp(argv[1], "--help") == 0)) )
	usage(argv);
    else if ( strcmp(argv[1], "--version") == 0 )
    {
	printf("%s %s\n", argv[0], VERSION);
	return EX_OK;
    }
    for (arg = 1; (arg < argc) && (*argv[arg] == '-'); ++arg)
    {
	if ( strcmp(argv[arg], "--verbose") == 0 )
	    fastq_trim_set_verbose(&tp, true);
	else if ( strcmp(argv[arg], "--exact-match") == 0 )
	    fastq_trim_set_adapter_match_function(&tp, bl_align_map_seq_exact);
	else if ( strcmp(argv[arg], "--3p-adapter1") == 0 )
	{
	    free(FASTQ_TRIM_ADAPTER1(&tp));
	    fastq_trim_set_adapter1(&tp, argv[++arg]);
	}
	else if ( strcmp(argv[arg], "--3p-adapter2") == 0 )
	{
	    free(FASTQ_TRIM_ADAPTER2(&tp));
	    fastq_trim_set_adapter2(&tp, argv[++arg]);
	}
	else if ( strcmp(argv[arg], "--min-match") == 0 )
	{
	    fastq_trim_set_min_match(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--max-mismatch-percent") == 0 )
	{
	    fastq_trim_set_max_mismatch_percent(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--min-qual") == 0 )
	{
	    fastq_trim_set_min_qual(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--min-length") == 0 )
	{
	    fastq_trim_set_min_length(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--phred-base") == 0 )
	{
	    fastq_trim_set_phred_base(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--polya-min-length") == 0 )
	{
	    fastq_trim_set_polya_min_len(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else
	    usage(argv);
    }
    xt_strupper(FASTQ_TRIM_ADAPTER1(&tp));
    
    fprintf(stdout, "\n*** FASTQ TRIM ***\n\n");
    fprintf(stdout, "  Minimum match:     %zu\n"
		    "  Minimum quality:   %u\n"
		    "  Minimum length:    %zu\n"
		    "  Phred base:        %u\n",
		    FASTQ_TRIM_MIN_MATCH(&tp),
		    FASTQ_TRIM_MIN_QUAL(&tp),
		    FASTQ_TRIM_MIN_LENGTH(&tp),
		    FASTQ_TRIM_PHRED_BASE(&tp));
    
    if ( FASTQ_TRIM_ADAPTER_MATCH_FUNCTION(&tp) == bl_align_map_seq_exact )
	fprintf(stdout, "  Adapter matching:  Exact\n");
    else
	fprintf(stdout, "  Adapter matching:  Smart\n"
			"  Maximum mismatch:  %u%%\n", FASTQ_TRIM_MAX_MISMATCH_PERCENT(&tp));
    if ( arg == argc )
	fprintf(stdout, "  Filename:          Standard input");
    else
	fprintf(stdout, "  Filename:          %s\n", argv[arg]);
    if ( arg + 2 < argc )
	fprintf(stdout, "  Filename:          %s\n", argv[arg+2]);

    if ( (status = fastq_trim_open_files(&tp, arg, argc, argv)) == EX_OK )
    {
	if ( FASTQ_TRIM_INFILE2(&tp) == NULL )
	    fastq_trim_single_reads(&tp);
	else
	    fastq_trim_paired_reads(&tp);
	fastq_trim_close_files(&tp);
	return EX_OK;
    }
    else if ( status == EX_USAGE )
	usage(argv);
    else
	return status;
}


void    usage(char *argv[])

{
    fprintf(stderr,
	    "Usage:\n\n"
	    "%s --help\n"
	    "%s --version\n"
	    "%s\n"
	    "   [--verbose]\n"
	    "   [--exact-match]\n"
	    "   [--3p-adapter1 seq]\n"
	    "   [--3p-adapter2 seq]\n"
	    "   [--min-match N]\n"
	    "   [--max-mismatch-percent N]\n"
	    "   [--min-qual N]\n"
	    "   [--min-length N]\n"
	    "   [--polya-min-length N]\n"
	    "   [--phred-base N]\n"
	    "   [infile1.fastq[.xz|.bz2|.gz]] [outfile1.fastq[.xz|.bz2|.gz]]\n\n"
	    "   [infile2.fastq[.xz|.bz2|.gz]] [outfile2.fastq[.xz|.bz2|.gz]]\n\n",
	    argv[0], argv[0], argv[0]);
    exit(EX_USAGE);
}
