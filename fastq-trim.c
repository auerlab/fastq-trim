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
size_t  bl_fastq_find_3p_adapter(const bl_fastq_t *read, const char *adapter, size_t min_match);
size_t  bl_fastq_find_3p_qual(const bl_fastq_t *read, unsigned min_qual, unsigned phred_base);
size_t  bl_fastq_3p_trim(bl_fastq_t *read, size_t new_len);

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
	else
	    usage(argv);
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
	index = bl_fastq_find_3p_adapter(&fastq_rec, adapter, min_match);
	if ( BL_FASTQ_SEQ_AE(&fastq_rec, index) != '\0' )
	{
	    ++adapters;
	    //fprintf(stderr, "%zu %s\n", BL_FASTQ_SEQ_LEN(&fastq_rec), BL_FASTQ_SEQ(&fastq_rec) + index);
	    bl_fastq_3p_trim(&fastq_rec, index);
	    //fprintf(stderr, "%zu %s\n", BL_FASTQ_SEQ_LEN(&fastq_rec), BL_FASTQ_SEQ(&fastq_rec) + index);
	    //getchar();
	}
	
	// FIXME: Support other PHRED bases
	index = bl_fastq_find_3p_qual(&fastq_rec, min_qual, 33);
	if ( BL_FASTQ_SEQ_AE(&fastq_rec, index) != '\0' )
	{
	    ++low_qual;
	    bl_fastq_3p_trim(&fastq_rec, index);
	}
	
	//fprintf(stderr, "%zu\n", BL_FASTQ_SEQ_LEN(&fastq_rec));
	if ( BL_FASTQ_SEQ_LEN(&fastq_rec) >= min_length )
	    bl_fastq_write(outstream, &fastq_rec, BL_FASTQ_SEQ_LEN(&fastq_rec));
	else
	    ++too_short;

	++records;
	if ( records % 100000 == 0 )
	{
	    fprintf(stderr,
		    "Reads: %lu  Adapters: %lu  Qual < %u: %lu  Len < %zu: %lu\r",
		    records, adapters, min_qual, low_qual, min_length, too_short);
	}
    }
    fprintf(stderr,
	    "Reads: %lu  Adapters: %lu  Qual < %u: %lu  Len < %zu: %lu\n",
	    records, adapters, min_qual, low_qual, min_length, too_short);

    bl_fastq_free(&fastq_rec);
    xt_fclose(outstream);
    xt_fclose(instream);
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
	    "    [infile.fastq[.xz|.bz2|.gz]] [outfile.fastq[.xz|.bz2|.gz]]\n\n",
	    argv[0], argv[0]);
    exit(EX_USAGE);
}

