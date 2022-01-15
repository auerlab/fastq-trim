
/*
 *  Generated by ./auto-gen-get-set
 *
 *  Mutator macros for setting with no sanity checking.  Use these to
 *  set structure members from functions outside the fastq_trim_t
 *  class.  These macros perform no data validation.  Hence, they achieve
 *  maximum performance where data are guaranteed correct by other means.
 *  Use the mutator functions (same name as the macro, but lower case)
 *  for more robust code with a small performance penalty.
 *
 *  These generated macros are not expected to be perfect.  Check and edit
 *  as needed before adding to your code.
 */

/* temp-fastq-trim-mutators.c */
int fastq_trim_set_verbose(fastq_trim_t *fastq_trim_ptr, _Bool new_verbose);
int fastq_trim_set_adapter_match_function(fastq_trim_t *fastq_trim_ptr, fastq_trim_afp_t new_adapter_match_function);
int fastq_trim_set_infile1(fastq_trim_t *fastq_trim_ptr, char *new_infile1);
int fastq_trim_set_infile1_ae(fastq_trim_t *fastq_trim_ptr, size_t c, char new_infile1_element);
int fastq_trim_set_infile1_cpy(fastq_trim_t *fastq_trim_ptr, char *new_infile1, size_t array_size);
int fastq_trim_set_outfile1(fastq_trim_t *fastq_trim_ptr, char *new_outfile1);
int fastq_trim_set_outfile1_ae(fastq_trim_t *fastq_trim_ptr, size_t c, char new_outfile1_element);
int fastq_trim_set_outfile1_cpy(fastq_trim_t *fastq_trim_ptr, char *new_outfile1, size_t array_size);
int fastq_trim_set_infile2(fastq_trim_t *fastq_trim_ptr, char *new_infile2);
int fastq_trim_set_infile2_ae(fastq_trim_t *fastq_trim_ptr, size_t c, char new_infile2_element);
int fastq_trim_set_infile2_cpy(fastq_trim_t *fastq_trim_ptr, char *new_infile2, size_t array_size);
int fastq_trim_set_outfile2(fastq_trim_t *fastq_trim_ptr, char *new_outfile2);
int fastq_trim_set_outfile2_ae(fastq_trim_t *fastq_trim_ptr, size_t c, char new_outfile2_element);
int fastq_trim_set_outfile2_cpy(fastq_trim_t *fastq_trim_ptr, char *new_outfile2, size_t array_size);
int fastq_trim_set_instream1(fastq_trim_t *fastq_trim_ptr, FILE *new_instream1);
int fastq_trim_set_instream1_ae(fastq_trim_t *fastq_trim_ptr, size_t c, FILE new_instream1_element);
int fastq_trim_set_instream1_cpy(fastq_trim_t *fastq_trim_ptr, FILE *new_instream1, size_t array_size);
int fastq_trim_set_outstream1(fastq_trim_t *fastq_trim_ptr, FILE *new_outstream1);
int fastq_trim_set_outstream1_ae(fastq_trim_t *fastq_trim_ptr, size_t c, FILE new_outstream1_element);
int fastq_trim_set_outstream1_cpy(fastq_trim_t *fastq_trim_ptr, FILE *new_outstream1, size_t array_size);
int fastq_trim_set_instream2(fastq_trim_t *fastq_trim_ptr, FILE *new_instream2);
int fastq_trim_set_instream2_ae(fastq_trim_t *fastq_trim_ptr, size_t c, FILE new_instream2_element);
int fastq_trim_set_instream2_cpy(fastq_trim_t *fastq_trim_ptr, FILE *new_instream2, size_t array_size);
int fastq_trim_set_outstream2(fastq_trim_t *fastq_trim_ptr, FILE *new_outstream2);
int fastq_trim_set_outstream2_ae(fastq_trim_t *fastq_trim_ptr, size_t c, FILE new_outstream2_element);
int fastq_trim_set_outstream2_cpy(fastq_trim_t *fastq_trim_ptr, FILE *new_outstream2, size_t array_size);
int fastq_trim_set_adapter1(fastq_trim_t *fastq_trim_ptr, char *new_adapter1);
int fastq_trim_set_adapter1_ae(fastq_trim_t *fastq_trim_ptr, size_t c, char new_adapter1_element);
int fastq_trim_set_adapter1_cpy(fastq_trim_t *fastq_trim_ptr, char *new_adapter1, size_t array_size);
int fastq_trim_set_adapter2(fastq_trim_t *fastq_trim_ptr, char *new_adapter2);
int fastq_trim_set_adapter2_ae(fastq_trim_t *fastq_trim_ptr, size_t c, char new_adapter2_element);
int fastq_trim_set_adapter2_cpy(fastq_trim_t *fastq_trim_ptr, char *new_adapter2, size_t array_size);
int fastq_trim_set_min_length(fastq_trim_t *fastq_trim_ptr, size_t new_min_length);
int fastq_trim_set_min_match(fastq_trim_t *fastq_trim_ptr, size_t new_min_match);
int fastq_trim_set_polya_min_len(fastq_trim_t *fastq_trim_ptr, size_t new_polya_min_len);
int fastq_trim_set_max_mismatch_percent(fastq_trim_t *fastq_trim_ptr, unsigned new_max_mismatch_percent);
int fastq_trim_set_min_qual(fastq_trim_t *fastq_trim_ptr, unsigned new_min_qual);
int fastq_trim_set_phred_base(fastq_trim_t *fastq_trim_ptr, unsigned new_phred_base);