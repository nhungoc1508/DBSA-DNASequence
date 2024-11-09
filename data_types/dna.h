#ifndef DNA_H
#define DNA_H

// #include "postgres.h"
// #include "fmgr.h"
#include "utils/builtins.h"
#include "varatt.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include "funcapi.h"
#include "executor/spi.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>

// #include "kmer.h"

#define VALID_NUCLEOTIDES   "ACGTacgt"
#define MAX_KMER_LEN        32

/******************************************************************************
 * TYPE STRUCT
 ******************************************************************************/

typedef struct {
    int32 vl_len_;
    int32 length;       
    char sequence[FLEXIBLE_ARRAY_MEMBER];  
} dna;

/******************************************************************************
 * AUXILIARY FUNCTIONS DECLARATION
 ******************************************************************************/

static dna *dna_make(int length, const char *sequence);
static bool is_valid_dna_sequence(const char *sequence);
static dna *dna_parse(const char *str);
static char *dna_to_string(const dna *seq);
static char *to_uppercase(char *data);

/******************************************************************************
 * AUXILIARY FUNCTIONS IMPLEMENTATION
 ******************************************************************************/

static dna *dna_make(int length, const char *sequence) {
    dna *result = (dna *) palloc(VARHDRSZ + sizeof(int32) + length + 1);
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + length + 1);
    result->length = length;
    // memcpy(result->sequence, sequence, length);   
    strcpy(result->sequence, sequence); 
    result->sequence[length] = '\0';  
    return result;
}

static bool is_valid_dna_sequence(const char *sequence) {
    for (int i = 0; sequence[i] != '\0'; i++) {
        if (strchr(VALID_NUCLEOTIDES, sequence[i]) == NULL) {
            return false;
        }
    }
    return true;
}

static dna *dna_parse(const char *str) {
    int length = strlen(str);

    if (!is_valid_dna_sequence(str)) {
        ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                        errmsg("Invalid character in DNA sequence")));
    }
    char *upper_str = to_uppercase(str);        
    return dna_make(length, upper_str);    
}

static char *dna_to_string(const dna *seq) {
    return psprintf("%s", seq->sequence);
}

#endif // DNA_H