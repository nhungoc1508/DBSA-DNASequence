// #include <math.h>
// #include <float.h>
// #include <stdlib.h>

// #include "postgres.h"
// #include "fmgr.h"
// #include "libpq/pqformat.h"
// #include "utils/fmgrprotos.h"

#include "dnasequence.h"
#include "spgist.h"

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

/**************kmer picksplit start******************/
static
Datum formKmerDatum(char *data, int k) {
    kmer *km = kmer_make(k, data);  /* Create kmer from data */
    return KmerPGetDatum(km);  /* Return it as a Datum */
}
  
 /*
  * Find the length of the common prefix of a and b
  */
 static int
 commonPrefix(const char *a, const char *b, int lena, int lenb)
 {
     int         i = 0;

     while (i < lena && i < lenb && nucleotide_matches(*a, *b))
     {  
         a++;
         b++;
         i++;
     }
  
     return i;
 }
 
/* Struct for sorting values in picksplit */
 typedef struct spgNodePtr
 {
     Datum       d;
     int         i;
     int16       c;
 } spgNodePtr;

/* qsort comparator to sort spgNodePtr structs by "c" */
 static int
 cmpNodePtr(const void *a, const void *b)
 {
     const spgNodePtr *aa = (const spgNodePtr *) a;
     const spgNodePtr *bb = (const spgNodePtr *) b;
  
     return pg_cmp_s16(aa->c, bb->c);
 }

Datum
 spg_kmer_picksplit(PG_FUNCTION_ARGS)
 {
     spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
     spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);
     kmer       *kmer0 = DatumGetKmerP(in->datums[0]);
     int         i,
                 commonLen;
     spgNodePtr *nodes;
  
     /* Identify longest common prefix, if any */
     commonLen = kmer0->k;
     for (i = 1; i < in->nTuples && commonLen > 0; i++)
     {
         kmer       *kmeri = DatumGetKmerP(in->datums[i]);
         int         tmp = commonPrefix(kmer0->data,
                                        kmeri->data,
                                        kmer0->k,
                                        kmeri->k);
  
         if (tmp < commonLen)
             commonLen = tmp;
     }
  
     /*
      * Limit the prefix length, if necessary, to ensure that the resulting
      * inner tuple will fit on a page.
      */
     /*commonLen = Min(commonLen, MAX_KMER_LEN);*/
  
     /* Set node prefix to be that string, if it's not empty */
     if (commonLen == 0)
     {
         out->hasPrefix = false;
     }
     else
     {
         out->hasPrefix = true;
         out->prefixDatum = formKmerDatum(kmer0->data, commonLen);
     }
  
     /* Extract the node label (first non-common byte) from each value */
     nodes = (spgNodePtr *) palloc(sizeof(spgNodePtr) * in->nTuples);
  
     for (i = 0; i < in->nTuples; i++)
     {
         kmer       *kmeri = DatumGetKmerP(in->datums[i]);
  
         if (commonLen < kmeri->k)
             nodes[i].c = *(unsigned char *) (kmer_i->data + commonLen);
         else
             nodes[i].c = -1;    /* use -1 if string is all common */
         nodes[i].i = i;
         nodes[i].d = in->datums[i];
     }
  
     /*
      * Sort by label values so that we can group the values into nodes.  This
      * also ensures that the nodes are ordered by label value, allowing the
      * use of binary search in searchChar.
      */
     qsort(nodes, in->nTuples, sizeof(*nodes), cmpNodePtr);
  
     /* And emit results */
     out->nNodes = 0;
     out->nodeLabels = (Datum *) palloc(sizeof(Datum) * in->nTuples);
     out->mapTuplesToNodes = (int *) palloc(sizeof(int) * in->nTuples);
     out->leafTupleDatums = (Datum *) palloc(sizeof(Datum) * in->nTuples);
  
     for (i = 0; i < in->nTuples; i++)
     {
         kmer       *kmeri = DatumGetKmerP(nodes[i].d);
         Datum       leafD;
  
         if (i == 0 || nodes[i].c != nodes[i - 1].c)
         {
             out->nodeLabels[out->nNodes] = Int16GetDatum(nodes[i].c);
             out->nNodes++;
         }
  
         if (commonLen < kmeri->k)
             leafD = formKmerDatum(kmeri->data + commonLen + 1,
                                   kmeri->k - commonLen - 1);
         else
             leafD = formKmerDatum(NULL, 0);
  
         out->leafTupleDatums[nodes[i].i] = leafD;
         out->mapTuplesToNodes[nodes[i].i] = out->nNodes - 1;
     }
  
     PG_RETURN_VOID();
 }
/**************kmer picksplit end******************/
 static inline int
 pg_cmp_s16(int16 a, int16 b)
 {
     return (int32) a - (int32) b;
 }
 
Datum
 spg_kmer_leaf_consistent(PG_FUNCTION_ARGS)
 {
     spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
     spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);
     int         level = in->level;
     text       *leafValue,
                *reconstrValue = NULL;
     char       *fullValue;
     int         fullLen;
     bool        res;
     int         j;
  
     /* all tests are exact */
     out->recheck = false;
  
     leafValue = DatumGetTextPP(in->leafDatum);
  
     /* As above, in->reconstructedValue isn't toasted or short. */
     if (DatumGetPointer(in->reconstructedValue))
         reconstrValue = (text *) DatumGetPointer(in->reconstructedValue);
  
     Assert(reconstrValue == NULL ? level == 0 :
            VARSIZE_ANY_EXHDR(reconstrValue) == level);
  
     /* Reconstruct the full string represented by this leaf tuple */
     fullLen = level + VARSIZE_ANY_EXHDR(leafValue);
     if (VARSIZE_ANY_EXHDR(leafValue) == 0 && level > 0)
     {
         fullValue = VARDATA(reconstrValue);
         out->leafValue = PointerGetDatum(reconstrValue);
     }
     else
     {
         text       *fullText = palloc(VARHDRSZ + fullLen);
  
         SET_VARSIZE(fullText, VARHDRSZ + fullLen);
         fullValue = VARDATA(fullText);
         if (level)
             memcpy(fullValue, VARDATA(reconstrValue), level);
         if (VARSIZE_ANY_EXHDR(leafValue) > 0)
             memcpy(fullValue + level, VARDATA_ANY(leafValue),
                    VARSIZE_ANY_EXHDR(leafValue));
         out->leafValue = PointerGetDatum(fullText);
     }
  
     /* Perform the required comparison(s) */
     res = true;
     for (j = 0; j < in->nkeys; j++)
     {
         StrategyNumber strategy = in->scankeys[j].sk_strategy;
         text       *query = DatumGetTextPP(in->scankeys[j].sk_argument);
         int         queryLen = VARSIZE_ANY_EXHDR(query);
         int         r;
  
         if (strategy == RTPrefixStrategyNumber)
         {
             /*
              * if level >= length of query then reconstrValue must begin with
              * query (prefix) string, so we don't need to check it again.
              */
             res = (level >= queryLen) ||
                 DatumGetBool(DirectFunctionCall2Coll(text_starts_with,
                                                      PG_GET_COLLATION(),
                                                      out->leafValue,
                                                      PointerGetDatum(query)));
  
             if (!res)           /* no need to consider remaining conditions */
                 break;
  
             continue;
         }
  
         if (SPG_IS_COLLATION_AWARE_STRATEGY(strategy))
         {
             /* Collation-aware comparison */
             strategy -= SPG_STRATEGY_ADDITION;
  
             /* If asserts enabled, verify encoding of reconstructed string */
             Assert(pg_verifymbstr(fullValue, fullLen, false));
  
             r = varstr_cmp(fullValue, fullLen,
                            VARDATA_ANY(query), queryLen,
                            PG_GET_COLLATION());
         }
         else
         {
             /* Non-collation-aware comparison */
             r = memcmp(fullValue, VARDATA_ANY(query), Min(queryLen, fullLen));
  
             if (r == 0)
             {
                 if (queryLen > fullLen)
                     r = -1;
                 else if (queryLen < fullLen)
                     r = 1;
             }
         }
  
         switch (strategy)
         {
             case BTLessStrategyNumber:
                 res = (r < 0);
                 break;
             case BTLessEqualStrategyNumber:
                 res = (r <= 0);
                 break;
             case BTEqualStrategyNumber:
                 res = (r == 0);
                 break;
             case BTGreaterEqualStrategyNumber:
                 res = (r >= 0);
                 break;
             case BTGreaterStrategyNumber:
                 res = (r > 0);
                 break;
             default:
                 elog(ERROR, "unrecognized strategy number: %d",
                      in->scankeys[j].sk_strategy);
                 res = false;
                 break;
         }
  
         if (!res)
             break;              /* no need to consider remaining conditions */
     }
  
     PG_RETURN_BOOL(res);
 }
