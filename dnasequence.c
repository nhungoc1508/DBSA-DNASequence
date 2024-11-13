// #include <math.h>
// #include <float.h>
// #include <stdlib.h>

// #include "postgres.h"
// #include "fmgr.h"
// #include "libpq/pqformat.h"
// #include "utils/fmgrprotos.h"

#include "dnasequence.h"

PG_MODULE_MAGIC;

/******************************************************************************
 * INPUT/OUTPUT ROUTINES
 ******************************************************************************/

// ********** dna **********
PG_FUNCTION_INFO_V1(dna_in);
PG_FUNCTION_INFO_V1(dna_out);

// ********** kmer **********
PG_FUNCTION_INFO_V1(kmer_in);
PG_FUNCTION_INFO_V1(kmer_out);
PG_FUNCTION_INFO_V1(kmer_cast_from_text);
PG_FUNCTION_INFO_V1(kmer_cast_to_text);

// ********** qkmer **********
PG_FUNCTION_INFO_V1(qkmer_in);
PG_FUNCTION_INFO_V1(qkmer_out);
PG_FUNCTION_INFO_V1(qkmer_cast_from_text);
PG_FUNCTION_INFO_V1(qkmer_cast_to_text);

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

// ********** dna **********
PG_FUNCTION_INFO_V1(dna_length);

// ********** kmer **********
PG_FUNCTION_INFO_V1(kmer_constructor);
PG_FUNCTION_INFO_V1(kmer_length);
PG_FUNCTION_INFO_V1(kmer_equals);
PG_FUNCTION_INFO_V1(kmer_starts_with);
PG_FUNCTION_INFO_V1(generate_kmers);
// Additional functions for the BTree operator class
PG_FUNCTION_INFO_V1(kmer_lt);
PG_FUNCTION_INFO_V1(kmer_le);
PG_FUNCTION_INFO_V1(kmer_gt);
PG_FUNCTION_INFO_V1(kmer_ge);
PG_FUNCTION_INFO_V1(kmer_cmp);

// ********** qkmer **********
PG_FUNCTION_INFO_V1(qkmer_constructor);
PG_FUNCTION_INFO_V1(qkmer_length);
PG_FUNCTION_INFO_V1(qkmer_contains);

/******************************************************************************
 * IMPLEMENTATION
 ******************************************************************************/

// ********** dna **********
Datum
dna_in(PG_FUNCTION_ARGS) {
    char *str = PG_GETARG_CSTRING(0);
    PG_RETURN_DNA_P(dna_parse(str));
}

Datum
dna_out(PG_FUNCTION_ARGS) {
    dna *seq = PG_GETARG_DNA_P(0);
    // char *result = dna_to_string(seq);
    char *result = strdup(seq->sequence);
    PG_FREE_IF_COPY(seq, 0);
    PG_RETURN_CSTRING(result);
}

Datum
dna_length(PG_FUNCTION_ARGS) {
    dna *seq = PG_GETARG_DNA_P(0);
    PG_RETURN_INT32(seq->length);
}

// ********** kmer **********
Datum
kmer_in(PG_FUNCTION_ARGS) {
    char *str = PG_GETARG_CSTRING(0);
    PG_RETURN_KMER_P(kmer_parse(&str));
}

Datum
kmer_out(PG_FUNCTION_ARGS) {
    kmer *c = PG_GETARG_KMER_P(0);
    char *result = strdup(c->data);
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_CSTRING(result);
}

Datum
kmer_cast_from_text(PG_FUNCTION_ARGS) {
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout, PointerGetDatum(txt)));
    PG_RETURN_KMER_P(kmer_parse(&str));
}

Datum
kmer_cast_to_text(PG_FUNCTION_ARGS) {
    kmer *c = PG_GETARG_KMER_P(0);
    text *out = (text *)DirectFunctionCall1(textin, PointerGetDatum(kmer_to_str(c)));
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_TEXT_P(out);
}

// Ref: https://doxygen.postgresql.org/varlena_8c.html#a7777d194920e57222a17b02166fc7232
Datum
kmer_constructor(PG_FUNCTION_ARGS) {
    int k = PG_GETARG_INT32(0);
    char *data = text_to_cstring(PG_GETARG_TEXT_PP(1));
    PG_RETURN_KMER_P(kmer_make(k, data));
}

Datum
kmer_length(PG_FUNCTION_ARGS) {
    kmer *c = PG_GETARG_KMER_P(0);
    int k = c->k;
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_INT32(k);
}

Datum
kmer_equals(PG_FUNCTION_ARGS) {
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    bool result;
    if (a->k != b->k) {
        result = false;
    } else {
        result = kmer_cmp_internal(a, b) == 0;
    }
    PG_FREE_IF_COPY(a, 0);
    PG_FREE_IF_COPY(b, 1);
    PG_RETURN_BOOL(result);
}

Datum
kmer_starts_with(PG_FUNCTION_ARGS) {
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

Datum
generate_kmers(PG_FUNCTION_ARGS) {    
    FuncCallContext *funcctx;
    int call_cntr;
    int max_calls;
    kmer *result_kmer;

    dna *dna_sequence = PG_GETARG_DNA_P(0);
    int k = PG_GETARG_INT32(1);
    
    if (SRF_IS_FIRSTCALL()) {        
        funcctx = SRF_FIRSTCALL_INIT();
        
        MemoryContext oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

        if (k <= 0 || k > dna_sequence->length || k > MAX_KMER_LEN) {
            ereport(ERROR, (errmsg("Invalid k value: must be between 1 and min(sequence length, %d)", MAX_KMER_LEN)));
        }

        funcctx->max_calls = dna_sequence->length - k + 1;
        funcctx->user_fctx = dna_sequence;

        MemoryContextSwitchTo(oldcontext);
    }

    funcctx = SRF_PERCALL_SETUP();

    call_cntr = funcctx->call_cntr;
    max_calls = funcctx->max_calls;

    if (call_cntr < max_calls) {        
        char *kmer_data = palloc(k + 1);
        strncpy(kmer_data, dna_sequence->sequence + call_cntr, k);
        kmer_data[k] = '\0';
        
        result_kmer = kmer_make(k, kmer_data);
        
        pfree(kmer_data);
        SRF_RETURN_NEXT(funcctx, KmerPGetDatum(result_kmer));
    } else {
        SRF_RETURN_DONE(funcctx);
    }
}

// BTree functions
Datum
kmer_lt(PG_FUNCTION_ARGS) {
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    bool result;
    if (a->k != b->k) {
        result = a->k < b->k;
    } else {
        result = kmer_cmp_internal(a, b) < 0;
    }
    PG_FREE_IF_COPY(a, 0);
    PG_FREE_IF_COPY(b, 1);
    PG_RETURN_BOOL(result);
}

Datum
kmer_le(PG_FUNCTION_ARGS) {
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    bool result;
    if (a->k != b->k) {
        result = a->k < b->k;
    } else {
        result = kmer_cmp_internal(a, b) <= 0;
    }
    PG_FREE_IF_COPY(a, 0);
    PG_FREE_IF_COPY(b, 1);
    PG_RETURN_BOOL(result);
}

Datum
kmer_gt(PG_FUNCTION_ARGS) {
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    bool result;
    if (a->k != b->k) {
        result = a->k > b->k;
    } else {
        result = kmer_cmp_internal(a, b) > 0;
    }
    PG_FREE_IF_COPY(a, 0);
    PG_FREE_IF_COPY(b, 1);
    PG_RETURN_BOOL(result);
}

Datum
kmer_ge(PG_FUNCTION_ARGS) {
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    bool result;
    if (a->k != b->k) {
        result = a->k > b->k;
    } else {
        result = kmer_cmp_internal(a, b) >= 0;
    }
    PG_FREE_IF_COPY(a, 0);
    PG_FREE_IF_COPY(b, 1);
    PG_RETURN_BOOL(result);
}

Datum
kmer_cmp(PG_FUNCTION_ARGS) {
    kmer *a = PG_GETARG_KMER_P(0);
    kmer *b = PG_GETARG_KMER_P(1);
    int result;
    if (a->k != b->k) {
        if (a->k > b->k) {
            result = 1;
        } else {
            result = -1;
        }
    } else {
        result = kmer_cmp_internal(a, b);
    }
    PG_FREE_IF_COPY(a, 0);
    PG_FREE_IF_COPY(b, 1);
    PG_RETURN_INT32(result);
}

// ********** qkmer **********
Datum
qkmer_in(PG_FUNCTION_ARGS) 
{
    char *str = PG_GETARG_CSTRING(0);
    PG_RETURN_QKMER_P(qkmer_parse(&str));
}

Datum
qkmer_out(PG_FUNCTION_ARGS) 
{
    qkmer *c = PG_GETARG_QKMER_P(0);
    char *result = strdup(c->data);
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_CSTRING(result);
}

Datum
qkmer_cast_from_text(PG_FUNCTION_ARGS) 
{
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout,
                 PointerGetDatum(txt)));
    PG_RETURN_QKMER_P(qkmer_parse(&str));
}

Datum 
qkmer_cast_to_text(PG_FUNCTION_ARGS) 
{
    qkmer *c = PG_GETARG_QKMER_P(0);
    text *out = (text *)DirectFunctionCall1(textin,
                 PointerGetDatum(qkmer_to_str(c)));
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_TEXT_P(out);
}

Datum
qkmer_constructor(PG_FUNCTION_ARGS) 
{
    int k = PG_GETARG_INT32(0);
    char *data = text_to_cstring(PG_GETARG_TEXT_PP(1));
    PG_RETURN_QKMER_P(qkmer_make(k, data));
}

Datum 
qkmer_length(PG_FUNCTION_ARGS) 
{
    qkmer *c = PG_GETARG_QKMER_P(0);
    int len = c->k;
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_INT32(len);
}

Datum
qkmer_contains(PG_FUNCTION_ARGS) {
    qkmer *pattern = PG_GETARG_QKMER_P(0);
    kmer *c = PG_GETARG_KMER_P(1);
    bool result;
    if (pattern->k != c->k) {
        result = false;
    } else {
        result = contains(pattern->data, c->data);
    }
    PG_FREE_IF_COPY(pattern, 0);
    PG_FREE_IF_COPY(c, 1);   
    PG_RETURN_BOOL(result);
}
