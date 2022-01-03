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
#include <xtend/file.h>
#include <biolibc/fastq.h>

#define ILLUMINA_UNIVERSAL      "AGATCGGAAGAG"
#define ILLUMINA_SMALL_RNA_3P   "TGGAATTCTCGG"
#define ILLUMINA_SMALL_RNA_5P   "GATCGTCGGACT"
#define NEXTERA                 "CTGTCTCTTATA"
#define SOLID                   "CGCCTTGGCCGT"
#define DEFAULT_ADAPTER         ILLUMINA_UNIVERSAL

void    usage(char *argv[]);
size_t  bl_find_3p_adapter(const bl_fastq_t *read, const char *adapter, size_t min_match);

int     main(int argc,char *argv[])

{
    unsigned long   records,
		    adapters,
		    too_short,
		    low_qual;
    char            *infile = NULL,
		    *outfile = NULL,
		    *adapter = DEFAULT_ADAPTER,
		    *end;
    size_t          index,
		    min_length = 30,
		    min_match = 3;
    unsigned        min_qual = 20;
    int             arg;
    FILE            *instream = stdin,
		    *outstream = stdout;
    bl_fastq_t      fastq_rec;

    if ( (argc == 2) && (strcmp(argv[1], "--help") == 0) )
	usage(argv);
    for (arg = 1; (arg < argc) && (*argv[arg] == '-'); ++arg)
    {
	if ( strcmp(argv[arg], "--3p-adapter") == 0 )
	{
	    adapter = argv[++arg];
	}
	else if ( strcmp(argv[arg], "--min-match") == 0 )
	{
	    min_match = strtoul(argv[++arg], &end, 10);
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--min-qual") == 0 )
	{
	    min_qual = strtoul(argv[++arg], &end, 10);
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[arg], "--min-length") == 0 )
	{
	    min_length = strtoul(argv[++arg], &end, 10);
	    if ( *end != '\0' )
		usage(argv);
	}
    }
    
    if ( arg < argc )
	infile = argv[arg++];
    if ( arg < argc )
	outfile = argv[arg];
    
    // Adapter used in our test file
    // FIXME: tolower() this ahead of time and tolower() the read bases
    // adapter="AGATCGGAAGAGCACAC";
    
    if ( infile != NULL )
	if ( (instream = xt_fopen(infile, "r")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], infile,
		    strerror(errno));
	    return EX_NOINPUT;
	}

    if ( outfile != NULL )
	if ( (outstream = xt_fopen(outfile, "w")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], infile,
		    strerror(errno));
	    return EX_CANTCREAT;
	}
    
    bl_fastq_init(&fastq_rec);
    records = adapters = too_short = low_qual = 0;
    while ( bl_fastq_read(instream, &fastq_rec) == BL_READ_OK )
    {
	index = bl_find_3p_adapter(&fastq_rec, adapter, min_match);
	if ( BL_FASTQ_SEQ_AE(&fastq_rec, index) != '\0' )
	{
	    ++adapters;
	    //bl_fastq_read_trim(&fastq_rec, index);
	}
	
	if ( BL_FASTQ_SEQ_LEN(&fastq_rec) >= min_length )
	    bl_fastq_write(outstream, &fastq_rec, BL_FASTQ_SEQ_LEN(&fastq_rec));
	else
	    ++too_short;

	++records;
	if ( records % 100000 == 0 )
	{
	    printf("Reads: %lu  Adapters: %lu\r", records, adapters);
	    fflush(stdout);
	}
    }
    printf("Reads: %lu  Adapters: %lu\n", records, adapters);

    bl_fastq_free(&fastq_rec);
    xt_fclose(outstream);
    xt_fclose(instream);
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s\n", argv[0]);
    exit(EX_USAGE);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/string.h>
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

size_t  bl_find_3p_adapter(const bl_fastq_t *read, const char *adapter,
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
