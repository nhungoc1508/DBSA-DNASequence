#include "postgres.h"
#include "fmgr.h"
#include "utils/builtins.h"
#include "varatt.h"

PG_MODULE_MAGIC;

#define VALID_NUCLEOTIDES "ACGTacgt"

typedef struct {
    int32 length;       
    char sequence[FLEXIBLE_ARRAY_MEMBER];  
} dna;

#define DatumGetDnaP(X) ((dna *) DatumGetPointer(X))
#define DnaPGetDatum(X) PointerGetDatum(X)
#define PG_GETARG_DNA_P(n) DatumGetDnaP(PG_GETARG_DATUM(n))
#define PG_RETURN_DNA_P(x) return DnaPGetDatum(x)

static dna *dna_make(int length, const char *sequence) {
    dna *result = (dna *) palloc(VARHDRSZ + length + 1);  
    SET_VARSIZE(result, VARHDRSZ + length + 1);  
    result->length = length;
    memcpy(result->sequence, sequence, length);  
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
    return dna_make(length, str);
}

static char *dna_to_string(const dna *seq) {
    return psprintf("%s", seq->sequence);
}


PG_FUNCTION_INFO_V1(dna_in);
Datum dna_in(PG_FUNCTION_ARGS) {
    char *str = PG_GETARG_CSTRING(0);
    PG_RETURN_DNA_P(dna_parse(str));
}

PG_FUNCTION_INFO_V1(dna_out);
Datum dna_out(PG_FUNCTION_ARGS) {
    dna *seq = PG_GETARG_DNA_P(0);
    char *result = dna_to_string(seq);
    PG_RETURN_CSTRING(result);
}

//length
PG_FUNCTION_INFO_V1(dna_length);
Datum dna_length(PG_FUNCTION_ARGS) {
    dna *seq = PG_GETARG_DNA_P(0);
    PG_RETURN_INT32(seq->length);
}
