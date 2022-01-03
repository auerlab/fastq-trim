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
#include <string.h>
#include <xtend/file.h>
#include <biolibc/fastq.h>

void    usage(char *argv[]);
char    *bl_find_adapter(const char *read, const char *adapter, size_t min_match);

int     main(int argc,char *argv[])

{
    unsigned long   records,
		    adapters;
    char    *infile,
	    *outfile,
	    *adapter,
	    *p;
    FILE    *instream,
	    *outstream;
    bl_fastq_t  fastq_rec;
    
    switch(argc)
    {
	case 4:
	    adapter = argv[1];
	    infile = argv[2];
	    outfile = argv[3];
	    break;
	
	default:
	    usage(argv);
    }

    // Adapter used in our test file
    // FIXME: tolower() this ahead of time and tolower() the read bases
    // adapter="AGATCGGAAGAGCACAC";
    
    if ( (instream = xt_fopen(infile, "r")) == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], infile,
		strerror(errno));
	return EX_NOINPUT;
    }
    
    if ( (outstream = xt_fopen(outfile, "w")) == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], infile,
		strerror(errno));
	return EX_CANTCREAT;
    }
    
    bl_fastq_init(&fastq_rec);
    records = adapters = 0;
    while ( bl_fastq_read(instream, &fastq_rec) == BL_READ_OK )
    {
	if ( (p = bl_find_adapter(BL_FASTQ_SEQ(&fastq_rec), adapter, 6)) != NULL )
	{
	    ++adapters;
	    *p = '\0';  // Trim adapter and everything after it
	}
	++records;
	bl_fastq_write(outstream, &fastq_rec, BL_FASTQ_SEQ_LEN(&fastq_rec));
	if ( records % 100000 == 0 )
	{
	    printf("Reads: %lu  Adapters: %lu\r", records, adapters);
	    fflush(stdout);
	}
    }

    printf("\nTotal reads processed = %lu\n", records);
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
 *      -lbiolibc
 *
 *  Description:
 *      Locate adapter (or a portion thereof if at the end) in string.
 *      Functions similarly to strstr(3).
 *  
 *  Arguments:
 *      read    Read sequence to be searched
 *      adapter Adapter sequence to be located
 *
 *  Returns:
 *      Pointer to located adapter sequence or NULL if not found
 *
 *  Examples:
 *
 *  Files:
 *
 *  Environment
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

char *bl_find_adapter(const char *read, const char *adapter, size_t min_match)

{
    const char *start, *pr, *pa;
    
    for (start = read; *start != '\0'; ++start)
    {
	for (pr = start, pa = adapter; *pr == *pa; ++pr, ++pa)
	    ;
	if ( ((*pr == '\0') && (pr - start >= min_match)) || (*pa == '\0') )
	    return (char *)start;
    }
    return NULL;
}
