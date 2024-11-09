#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "varatt.h"

PG_MODULE_MAGIC;

#define VALID_NUCLEOTIDES "ACGTacgt"

typedef struct {
    int32 vl_len_;        // header          
    int32 length;         // length of dna        
    char data[FLEXIBLE_ARRAY_MEMBER]; 
} dna;

Datum dna_in(PG_FUNCTION_ARGS);
Datum dna_out(PG_FUNCTION_ARGS);
Datum dna_length(PG_FUNCTION_ARGS);
static dna *dna_make(const char *input);

PG_FUNCTION_INFO_V1(dna_in);
PG_FUNCTION_INFO_V1(dna_out);
PG_FUNCTION_INFO_V1(dna_length);

static bool is_valid_dna_sequence(const char *sequence) {
    for (int i = 0; sequence[i] != '\0'; i++) {
        if (strchr(VALID_NUCLEOTIDES, sequence[i]) == NULL) {
            return false;
        }
    }
    return true;
}

static dna *dna_make(const char *input) {
    if (!is_valid_dna_sequence(input))
        ereport(ERROR, (errmsg("VALID_NUCLEOTIDES")));

    int len = strlen(input);
        
    dna *result = (dna *) palloc(VARHDRSZ + sizeof(int32) + len); // VARHDRSZ for 'vl_len_'
    SET_VARSIZE(result, VARHDRSZ + sizeof(int32) + len); // Total size includes length and data fields

    result->length = len; 
    memcpy(result->data, input, len); 

    return result;
}

Datum
dna_in(PG_FUNCTION_ARGS) {
    char *input = PG_GETARG_CSTRING(0);
    dna *result = dna_make(input);
    PG_RETURN_POINTER(result);
}

Datum
dna_out(PG_FUNCTION_ARGS) {
    dna *dna_seq = (dna *) PG_GETARG_POINTER(0);
    char *output = (char *) palloc(dna_seq->length + 1); 
    memcpy(output, dna_seq->data, dna_seq->length); 
    output[dna_seq->length] = '\0'; // Null-terminate the string

    PG_RETURN_CSTRING(output);
}

Datum
dna_length(PG_FUNCTION_ARGS) {
    dna *dna_seq = (dna *) PG_GETARG_POINTER(0);
    PG_RETURN_INT32(dna_seq->length); 
}
