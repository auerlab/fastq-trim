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
#include "trim.h"

void    usage(char *argv[]);

int     main(int argc,char *argv[])

{
    int     arg, status;
    char    *end;
    trim_t  tp;
    
    trim_init(&tp);
    
    if ( (argc == 2) && (strcmp(argv[1], "--help") == 0) )
	usage(argv);
    for (arg = 1; (arg < argc) && (*argv[arg] == '-'); ++arg)
    {
	if ( strcmp(argv[arg], "--verbose") == 0 )
	    trim_set_verbose(&tp, true);
	if ( strcmp(argv[arg], "--exact-match") == 0 )
	    trim_set_adapter_match_function(&tp, bl_fastq_find_adapter_exact);
	else if ( strcmp(argv[arg], "--3p-adapter1") == 0 )
	{
	    free(TRIM_ADAPTER1(&tp));
	    trim_set_adapter1(&tp, argv[++arg]);
	}
	else if ( strcmp(argv[arg], "--3p-adapter2") == 0 )
	{
	    free(TRIM_ADAPTER2(&tp));
	    trim_set_adapter2(&tp, argv[++arg]);
	}
	else if ( strcmp(argv[arg], "--min-match") == 0 )
	{
	    trim_set_min_match(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--max-mismatch-percent") == 0 )
	{
	    trim_set_max_mismatch_percent(&tp, strtoul(argv[++arg], &end, 10));
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
	else if ( strcmp(argv[arg], "--phred-base") == 0 )
	{
	    trim_set_phred_base(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--polya-min-length") == 0 )
	{
	    trim_set_polya_min_len(&tp, strtoul(argv[++arg], &end, 10));
	    if ( *end != '\0' )
		usage(argv);
	}
	else
	    usage(argv);
    }
    strupper(TRIM_ADAPTER1(&tp));
    
    fprintf(stderr, "\n*** FASTQ TRIM ***\n\n");
    fprintf(stderr, "  Minimum match:     %zu\n"
		    "  Minimum quality:   %u\n"
		    "  Minimum length:    %zu\n"
		    "  Phred base:        %u\n",
		    TRIM_MIN_MATCH(&tp),
		    TRIM_MIN_QUAL(&tp),
		    TRIM_MIN_LENGTH(&tp),
		    TRIM_PHRED_BASE(&tp));
    
    if ( TRIM_ADAPTER_MATCH_FUNCTION(&tp) == bl_fastq_find_adapter_exact )
	fprintf(stderr, "  Adapter matching:  Exact\n");
    else
	fprintf(stderr, "  Adapter matching:  Smart\n"
			"  Maximum mismatch:  %u%%\n", TRIM_MAX_MISMATCH_PERCENT(&tp));
    if ( arg == argc )
	fprintf(stderr, "  Filename:          Standard input");
    else
	fprintf(stderr, "  Filename:          %s\n", argv[arg]);
    if ( arg + 2 < argc )
	fprintf(stderr, "  Filename:          %s\n", argv[arg+2]);

    if ( (status = trim_open_files(&tp, arg, argc, argv)) == EX_OK )
    {
	if ( TRIM_INFILE2(&tp) == NULL )
	    trim_single_reads(&tp);
	else
	    trim_paired_reads(&tp);
	trim_close_files(&tp);
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
	    "%s\n"
	    "   [--verbose]\n"
	    "   [--exact-match]\n"
	    "   [--3p-adapter1 seq]\n"
	    "   [--3p-adapter2 seq]\n"
	    "   [--min-match N]\n"
	    "   [--max-mismatch-percent N]\n"
	    "   [--min-qual N]\n"
	    "   [--min-length N]\n"
	    "   [--polya_min-length N]\n"
	    "   [--phred-base N]\n"
	    "   [infile1.fastq[.xz|.bz2|.gz]] [outfile1.fastq[.xz|.bz2|.gz]]\n\n"
	    "   [infile2.fastq[.xz|.bz2|.gz]] [outfile2.fastq[.xz|.bz2|.gz]]\n\n",
	    argv[0], argv[0]);
    exit(EX_USAGE);
}
