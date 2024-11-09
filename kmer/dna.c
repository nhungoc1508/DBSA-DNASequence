#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"

PG_MODULE_MAGIC;

#define MAX_KMER_LEN        32
#define VALID_NUCLEOTIDES   4
char nucleotides[VALID_NUCLEOTIDES * 2 + 1] = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', '\0'};

/*****************************************************************************/

typedef struct {
    int k;
    char *data;
} kmer;

#define DatumGetKmerP(X) ((kmer *) DatumGetPointer(X))
#define KmerPGetDatum(X) PointerGetDatum(X)
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x)

/*****************************************************************************/

static kmer *kmer_make(int k, char *data) {
    kmer *c = palloc0(sizeof(kmer));
    c->k = k;
    c->data = strdup(data);
    return c;
}

/*****************************************************************************/

/* 
 * Remove leading whitespaces etc
 */
static void p_whitespace(char **str) {
  while (**str == ' ' || **str == '\n' || **str == '\r' || **str == '\t')
    *str += 1;
}

static void ensure_end_input(char **str, bool end) {
    if (end) {
        p_whitespace(str);
        if (**str != 0) {
            ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION), errmsg("Could not parse temporal value")));
        }
    }
}

static bool is_valid_nucleotide(char *str) {
    for (int i=0; i<8; i++) {
        if (*str == nucleotides[i]) {
            return true;
        }
    }
    return false;
}

static bool is_valid_kmer(char **str) {
    char *alt_str = *str;
    while (*alt_str) {
        if (is_valid_nucleotide(alt_str)) {
            alt_str += 1;
        } else {
            return false;
        }
    }
    return true;
}

static char *to_uppercase(char *data) {
    int length = strlen(data);
    char *upper_str = malloc(length + 1);
    for (int i=0; data[i] != '\0'; i++) {
        if (data[i] >= 'a' && data[i] <= 'z') {
            upper_str[i] = data[i] - 32;
        } else {
            upper_str[i] = data[i];
        }
    }
    upper_str[length] = '\0';
    return upper_str;
}

static kmer *kmer_parse(char **str) {
    int k;
    char *data;
    bool ret = is_valid_kmer(str);
    if (!ret)
        ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION), errmsg("Invalid input syntax for type kmer")));
    // data = strdup(*str);
    char *upper_str;
    upper_str = to_uppercase(*str);
    data = strdup(upper_str);
    k = (int)strlen(data);
    if (k > MAX_KMER_LEN)
        ereport(ERROR, (errcode(ERRCODE_STRING_DATA_RIGHT_TRUNCATION), errmsg("Input exceeds maximum length allowed for type kmer (32)")));
    return kmer_make(k, data);
}

static char *kmer_to_str(const kmer *c) {
    char *result = strdup(c->data);
    return result;
}

static bool starts_with(char *prefix, char *c) {
    bool result;
    result = strncmp(prefix, c, strlen(prefix)) == 0;
    return result;
}

/*****************************************************************************/

PG_FUNCTION_INFO_V1(kmer_in);
Datum kmer_in(PG_FUNCTION_ARGS) {
    char *str = PG_GETARG_CSTRING(0);
    PG_RETURN_KMER_P(kmer_parse(&str));
}

PG_FUNCTION_INFO_V1(kmer_out);
Datum kmer_out(PG_FUNCTION_ARGS) {
    kmer *c = PG_GETARG_KMER_P(0);
    char *result = strdup(c->data);
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_CSTRING(result);
}

// PG_FUNCTION_INFO_V1(kmer_recv);
// Datum kmer_recv(PG_FUNCTION_ARGS) {
//     // Todo
// }

// PG_FUNCTION_INFO_V1(kmer_send);
// Datum kmer_send(PG_FUNCTION_ARGS) {
//     // Todo
// }

PG_FUNCTION_INFO_V1(kmer_cast_from_text);
Datum kmer_cast_from_text(PG_FUNCTION_ARGS) {
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout, PointerGetDatum(txt)));
    PG_RETURN_KMER_P(kmer_parse(&str));
}

PG_FUNCTION_INFO_V1(kmer_cast_to_text);
Datum kmer_cast_to_text(PG_FUNCTION_ARGS) {
    kmer *c = PG_GETARG_KMER_P(0);
    text *out = (text *)DirectFunctionCall1(textin, PointerGetDatum(kmer_to_str(c)));
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_TEXT_P(out);
}

/******************************************************************************
 * CONSTRUCTOR
 ******************************************************************************/
// Ref: https://doxygen.postgresql.org/varlena_8c.html#a7777d194920e57222a17b02166fc7232
PG_FUNCTION_INFO_V1(kmer_constructor);
Datum kmer_constructor(PG_FUNCTION_ARGS) {
    int k = PG_GETARG_INT32(0);
    char *data = text_to_cstring(PG_GETARG_TEXT_PP(1));
    PG_RETURN_KMER_P(kmer_make(k, data));
}

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

PG_FUNCTION_INFO_V1(kmer_length);
Datum kmer_length(PG_FUNCTION_ARGS) {
    kmer *c = PG_GETARG_KMER_P(0);
    int k = c->k;
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_INT32(k);
}

PG_FUNCTION_INFO_V1(kmer_equals);
Datum kmer_equals(PG_FUNCTION_ARGS) {
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    bool result;
    if (a->k != b->k) {
        result = false;
    } else {
        result = strcmp(a->data, b->data) == 0;
    }
    PG_FREE_IF_COPY(a, 0);
    PG_FREE_IF_COPY(b, 1);
    PG_RETURN_BOOL(result);
}

PG_FUNCTION_INFO_V1(kmer_starts_with);
Datum kmer_starts_with(PG_FUNCTION_ARGS) {
    kmer *prefix = PG_GETARG_KMER_P(0);
    kmer *c = PG_GETARG_KMER_P(1);
    bool result;
    if (prefix->k > c->k) {
        result = false;
    } else {
        result = starts_with(prefix->data, c->data);
    }
    PG_FREE_IF_COPY(prefix, 0);
    PG_FREE_IF_COPY(c, 1);   
    PG_RETURN_BOOL(result);
}