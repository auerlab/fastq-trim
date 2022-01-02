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

int     main(int argc,char *argv[])

{
    unsigned long   records,
		    adapters;
    char    *infile,
	    *adapter;
    FILE    *instream;
    bl_fastq_t  fastq_rec;
    
    switch(argc)
    {
	case 2:
	    infile = argv[1];
	    break;
	
	default:
	    usage(argv);
    }
    
    adapter="AGATCGGAAGAGCACAC";
    if ( (instream = xt_fopen(infile, "r")) == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], infile,
		strerror(errno));
	return EX_NOINPUT;
    }
    
    bl_fastq_init(&fastq_rec);
    records = adapters = 0;
    while ( bl_fastq_read(instream, &fastq_rec) == BL_READ_OK )
    {
	if ( strstr(BL_FASTQ_SEQ(&fastq_rec), adapter) != NULL )
	    ++adapters;
	++records;
	if ( records % 100000 == 0 )
	{
	    printf("Reads: %lu  Adapters: %lu\r", records, adapters);
	    fflush(stdout);
	}
    }

    printf("\nTotal reads processed = %lu\n", records);
    bl_fastq_free(&fastq_rec);
    xt_fclose(instream);
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s\n", argv[0]);
    exit(EX_USAGE);
}
