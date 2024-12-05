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
PG_FUNCTION_INFO_V1(kmer_starts_with_swapped);
PG_FUNCTION_INFO_V1(generate_kmers);
// Additional functions for the BTree operator class
PG_FUNCTION_INFO_V1(kmer_lt);
PG_FUNCTION_INFO_V1(kmer_le);
PG_FUNCTION_INFO_V1(kmer_gt);
PG_FUNCTION_INFO_V1(kmer_ge);
PG_FUNCTION_INFO_V1(kmer_cmp);
// Additional functions for the hash operator class
PG_FUNCTION_INFO_V1(kmer_hash);
// Additional functions for the SP-GiST operator class
PG_FUNCTION_INFO_V1(spgist_kmer_config);
PG_FUNCTION_INFO_V1(spgist_kmer_choose);
PG_FUNCTION_INFO_V1(spgist_kmer_picksplit);
PG_FUNCTION_INFO_V1(spgist_kmer_inner_consistent);
PG_FUNCTION_INFO_V1(spgist_kmer_leaf_consistent);

// ********** qkmer **********
PG_FUNCTION_INFO_V1(qkmer_constructor);
PG_FUNCTION_INFO_V1(qkmer_length);
PG_FUNCTION_INFO_V1(qkmer_contains);
PG_FUNCTION_INFO_V1(qkmer_contains_swapped);

/******************************************************************************
 * IMPLEMENTATION
 ******************************************************************************/

// ********** dna **********
Datum
dna_in(PG_FUNCTION_ARGS) {
    char *str = PG_GETARG_CSTRING(0);
    if (str == NULL || strlen(str) == 0) {
        ereport(ERROR,
                (errcode(ERRCODE_NULL_VALUE_NOT_ALLOWED),
                 errmsg("Input to dna type cannot be NULL or empty")));
    }
    PG_RETURN_DNA_P(dna_parse(str));
}

Datum
dna_out(PG_FUNCTION_ARGS) {
    dna *seq = PG_GETARG_DNA_P(0);
    char *result = pstrdup(seq->sequence);
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
    char *result = pstrdup(c->data);
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
kmer_starts_with_swapped(PG_FUNCTION_ARGS) {
    kmer *c = PG_GETARG_KMER_P(0);
    kmer *prefix = PG_GETARG_KMER_P(1);
    bool result;
    if (prefix->k > c->k) {
        result = false;
    } else {
        result = starts_with(prefix->data, c->data);
    }
    PG_FREE_IF_COPY(prefix, 1);
    PG_FREE_IF_COPY(c, 0);   
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
        funcctx->user_fctx = palloc(sizeof(dna));
        memcpy(funcctx->user_fctx, dna_sequence, sizeof(dna));

        MemoryContextSwitchTo(oldcontext);
    }

    funcctx = SRF_PERCALL_SETUP();

    call_cntr = funcctx->call_cntr;
    max_calls = funcctx->max_calls;

    if (call_cntr < max_calls) {        
        dna *sequence = (dna *)funcctx->user_fctx;
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

// Hash function
Datum
kmer_hash(PG_FUNCTION_ARGS) {
    kmer *a = PG_GETARG_KMER_P(0);
    int32 hash;
    hash = (int32) hash_any((unsigned char *) a->data, strlen(a->data));
    PG_FREE_IF_COPY(a, 0);
    PG_RETURN_INT32(hash);
}
// SP-GiST functions
Datum
spgist_kmer_config(PG_FUNCTION_ARGS) {
    // spgConfigIn *cfgin = (spgConfigIn *) PG_GETARG_POINTER(0);
    spgConfigOut *cfg = (spgConfigOut *) PG_GETARG_POINTER(1);

    // cfgin->attType = TEXTOID; //type of data index will store (kmer is text)
    cfg->prefixType = TEXTOID;
    cfg->labelType = INT2OID; // labels determine how data is partitioned
    // cfg->leafType = TEXTOID;
    cfg->canReturnData = true; // true so index can return data when queried
    cfg->longValuesOK = false; 
    PG_RETURN_VOID();
}

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
    int            i = 0;

    while (i < lena && i < lenb && nucleotide_matches(*a, *b))
    {
        a++;
        b++;
        i++;
    }

    return i;
}

/*
 * Binary search an array of int16 datums for a match to c
 *
 * On success, *i gets the match location; on failure, it gets where to insert
 */
static bool
searchChar(Datum *nodeLabels, int nNodes, int16 c, int *i)
{
    int            StopLow = 0,
                StopHigh = nNodes;

    while (StopLow < StopHigh)
    {
        int            StopMiddle = (StopLow + StopHigh) >> 1;
        int16        middle = DatumGetInt16(nodeLabels[StopMiddle]);

        if (c < middle)
            StopHigh = StopMiddle;
        else if (c > middle)
            StopLow = StopMiddle + 1;
        else
        {
            *i = StopMiddle;
            return true;
        }
    }

    *i = StopHigh;
    return false;
}

Datum
spgist_kmer_choose(PG_FUNCTION_ARGS) {
    spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);

    kmer    *inKmer = DatumGetKmerP(in->datum);
    char    *inData = inKmer->data;
    int      inK = inKmer->k;
    char    *prefixData = NULL;
    int      prefixK = 0;
    int      commonLen = 0;
    int16    nodeChar = 0;
    int      i = 0;

    /* Check for prefix match, set nodeChar to first byte after prefix */
    if (in->hasPrefix)
    {
        kmer *prefixKmer = DatumGetKmerP(in->prefixDatum);

        prefixData = prefixKmer->data;
        prefixK = prefixKmer->k;

        commonLen = commonPrefix(inData + in->level,
                                 prefixData,
                                 inK - in->level,
                                 prefixK);

        if (commonLen == prefixK)
        {
            /* node label --- first non-common character */
            if (inK - in->level > commonLen)
                nodeChar = *(unsigned char *) (inData + in->level + commonLen);
            else
                nodeChar = -1; /* completely common values */
        }
        else /* Not match prefix  -> split tuple */
        {
            out->resultType = spgSplitTuple;

            if (commonLen == 0)
            {
                out->result.splitTuple.prefixHasPrefix = false;
            }
            else
            {
                out->result.splitTuple.prefixHasPrefix = true;
                out->result.splitTuple.prefixPrefixDatum =
                    formKmerDatum(prefixData, commonLen);
            }
            out->result.splitTuple.prefixNNodes = 1;
            out->result.splitTuple.prefixNodeLabels =
                (Datum *) palloc(sizeof(Datum));
            out->result.splitTuple.prefixNodeLabels[0] =
                Int16GetDatum(*(unsigned char *) (prefixData + commonLen));

            out->result.splitTuple.childNodeN = 0;

            if (prefixK - commonLen == 1)
            {
                out->result.splitTuple.postfixHasPrefix = false;
            }
            else
            {
                out->result.splitTuple.postfixHasPrefix = true;
                out->result.splitTuple.postfixPrefixDatum =
                    formKmerDatum(prefixData + commonLen + 1, prefixK - commonLen - 1);
            }

            PG_RETURN_VOID();
        }
    }
    else if (inK > in->level)
    {
        nodeChar = *(unsigned char *) (inData + in->level); /* node label = 1st character after the current level */
    }
    else
    {
        nodeChar = -1;
    }

    /* Look up nodeChar in the node label array */
    if (searchChar(in->nodeLabels, in->nNodes, nodeChar, &i))
    {
        /* Descend to the existing node */
        int levelAdd;

        out->resultType = spgMatchNode;
        out->result.matchNode.nodeN = i;
        levelAdd = commonLen;
        if (nodeChar >= 0)
            levelAdd++;
        out->result.matchNode.levelAdd = levelAdd;
        if (inK - in->level - levelAdd > 0)
        {
            out->result.matchNode.restDatum =
                formKmerDatum(inData + in->level + levelAdd,
                              inK - in->level - levelAdd);
        }
        else
        {
            out->result.matchNode.restDatum = formKmerDatum(NULL, 0);
        }
    }
    else if (in->allTheSame)
    {
        /* Cannot use AddNode; split the tuple */
        out->resultType = spgSplitTuple;
        out->result.splitTuple.prefixHasPrefix = in->hasPrefix;
        out->result.splitTuple.prefixPrefixDatum = in->prefixDatum;
        out->result.splitTuple.prefixNNodes = 1;
        out->result.splitTuple.prefixNodeLabels = (Datum *) palloc(sizeof(Datum));
        out->result.splitTuple.prefixNodeLabels[0] = Int16GetDatum(-2);
        out->result.splitTuple.childNodeN = 0;
        out->result.splitTuple.postfixHasPrefix = false;
    }
    else
    {
        /* Add a new node for the unseen nodeChar */
        out->resultType = spgAddNode;
        out->result.addNode.nodeLabel = Int16GetDatum(nodeChar);
        out->result.addNode.nodeN = i;
    }

    PG_RETURN_VOID();
}

static inline int
pg_cmp_s16(int16 a, int16 b)
{
   return (int32) a - (int32) b;
}

/* qsort comparator to sort spgNodePtr structs by "c" */
static int
cmpNodePtr(const void *a, const void *b)
{
    const spgNodePtr *aa = (const spgNodePtr *) a;
    const spgNodePtr *bb = (const spgNodePtr *) b;

    return pg_cmp_s16(aa->c, bb->c);
}

Datum
spgist_kmer_picksplit(PG_FUNCTION_ARGS) {
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
            nodes[i].c = *(unsigned char *) (kmeri->data + commonLen);
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

Datum
spgist_kmer_inner_consistent(PG_FUNCTION_ARGS){
    // 
    spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
    // elog(NOTICE, "inner called, nkeys: %d", in->nkeys);
    spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);
    // bool        collate_is_c = pg_newlocale_from_collation(PG_GET_COLLATION())->collate_is_c;
    kmer       *reconstructedValue;
    kmer       *reconstrKmer;
    int         maxReconstrLen;
    kmer       *prefixKmer = NULL;
    int         prefixSize = 0;
    int         i;

    /*
    * Reconstruct values represented at this tuple, including parent data,
    * prefix of this tuple if any, and the node label if it's non-dummy.
    * in->level should be the length of the previously reconstructed value,
    * and the number of bytes added here is prefixSize or prefixSize + 1.
    *
    * Note: we assume that in->reconstructedValue isn't toasted and doesn't
    * have a short varlena header.  This is okay because it must have been
    * created by a previous invocation of this routine, and we always emit
    * long-format reconstructed values.
    */
    reconstructedValue = (kmer *) DatumGetPointer(in->reconstructedValue);

    // Assert(reconstructedValue == NULL ? in->level == 0 :
    //       (reconstructedValue->k) == in->level);

    maxReconstrLen = in->level + 1;

    if (in->hasPrefix) {
        prefixKmer = (kmer *) DatumGetPointer(in->prefixDatum);
        prefixSize = (prefixKmer->k);
        maxReconstrLen += prefixSize;
    }
    // reconstrKmer = new reconstruction buffer
    // The previously reconstructed value and any prefix
    // will be copied to the new reconstruction buffer.
    // elog(NOTICE, "[inner] reconstrKmer maxReconstrLen: %d, palloc: %d, pfsize: %d, in->level: %d", maxReconstrLen, VARHDRSZ + sizeof(int32) + maxReconstrLen + 1, prefixSize, in->level);
    reconstrKmer = (kmer *) palloc(VARHDRSZ + sizeof(int32) + maxReconstrLen + 1);
    SET_VARSIZE(reconstrKmer, VARHDRSZ + sizeof(int32) + maxReconstrLen + 1);
    reconstrKmer->k = maxReconstrLen;

    if (in->level)
        memcpy((reconstrKmer->data),
               (reconstructedValue->data),
               in->level);
    if (prefixSize)
        memcpy(((char *) (reconstrKmer->data)) + in->level,
                (prefixKmer->data),
                prefixSize);
    /* last byte of reconstrKmer will be filled in below */
    /*
     * Scan the child nodes.  For each one, complete the reconstructed value
     * and see if it's consistent with the query.  If so, emit an entry into
     * the output arrays.
     */
    out->nodeNumbers = (int *) palloc(sizeof(int) * in->nNodes);
    out->levelAdds = (int *) palloc(sizeof(int) * in->nNodes);
    out->reconstructedValues = (Datum *) palloc(sizeof(Datum) * in->nNodes);
    out->nNodes = 0;

    for (i = 0; i < in->nNodes; i++) {
        int16   nodeChar = DatumGetInt16(in->nodeLabels[i]);
        int     thisLen;
        bool    res = true;
        int     j;

        /* If nodeChar is a dummy value, don't include it in data */
        if (nodeChar <= 0) {
            thisLen = maxReconstrLen - 1;
        } else {
            ((unsigned char *) (reconstrKmer->data))[maxReconstrLen - 1] = nodeChar;
            thisLen = maxReconstrLen;
        }
        ((unsigned char *) (reconstrKmer->data))[thisLen] = '\0';
        for (j = 0; j < in->nkeys; j++) {
            StrategyNumber strategy = in->scankeys[j].sk_strategy;
            kmer    *inKmer;
            qkmer   *inQkmer;
            int      inSize;
            int      r;

            /*
            * If it's a collation-aware operator, but the collation is C, we
            * can treat it as non-collation-aware.  With non-C collation we
            * need to traverse whole tree :-( so there's no point in making
            * any check here.  (Note also that our reconstructed value may
            * well end with a partial multibyte character, so that applying
            * any encoding-sensitive test to it would be risky anyhow.)
            */

            inKmer = DatumGetKmerP(in->scankeys[j].sk_argument);
            inSize = (inKmer->k);

            size_t minLen = Min(inSize, thisLen);

            
            switch (strategy)
            {
                case BTEqualStrategyNumber:
                    r = 0;
                    for (size_t i = 0; i < minLen; i++) {
                        if (!nucleotide_matches(reconstrKmer->data[i], inKmer->data[i])) {
                            r = reconstrKmer->data[i] - inKmer->data[i];
                            break;
                        }
                    }
                    if (r != 0 || inSize < thisLen)
                        res = false;
                    break;
                case BTStartsWithStrategyNumber:
                    r = 0;
                    r = strncmp(inKmer->data, reconstrKmer->data, minLen);
                    res = (r == 0);
                    break;
                case BTContainsStrategyNumber:
                    inQkmer = DatumGetQkmerP(in->scankeys[j].sk_argument);
                    inSize = inQkmer->k;
                    minLen = Min(inSize, thisLen);
                    r = 0;
                    for (size_t i = 0; i < minLen; i++) {
                        if (!nucleotide_matches(inQkmer->data[i], reconstrKmer->data[i])) {
                            r = reconstrKmer->data[i] - inKmer->data[i];
                            break;
                        }
                    }
                    if (r != 0 || inSize < thisLen)
                        res = false;
                    // elog(NOTICE, "[inner] CONTAINS: inQkmer: %s, k: %d, minLen: %d, reconst: %s, inSize: %d, thisLen: %d, r: %d", inQkmer->data, inQkmer->k, minLen, reconstrKmer->data, inSize, thisLen, r);
                    break;
            }

            if (!res)
                break;          /* no need to consider remaining conditions */
        }

        if (res)
        {
            out->nodeNumbers[out->nNodes] = i;
            out->levelAdds[out->nNodes] = thisLen - in->level;
            SET_VARSIZE(reconstrKmer, VARHDRSZ + sizeof(int32) + thisLen + 1);
            // elog(NOTICE, "[inner] out node: %s", reconstrKmer->data);
            // if (thisLen == 32) {
            //     elog(NOTICE, "[inner], reconstrKmer: %s", reconstrKmer->data);
            // }
            out->reconstructedValues[out->nNodes] =
                datumCopy(KmerPGetDatum(reconstrKmer), false, -1);
            out->nNodes++;
        }
    }
    PG_RETURN_VOID();
}

Datum
spgist_kmer_leaf_consistent(PG_FUNCTION_ARGS)
{
    spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
    // elog(NOTICE, "leaf called, nkeys: %d", in->nkeys);
    spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);
    int         level = in->level;
    kmer       *leafValue,
               *reconstrValue = NULL;
    char       *fullValue;
    int         fullLen;
    bool        res;
    int         j;

    out->recheck = false;

    leafValue = DatumGetKmerP(in->leafDatum);

    /* As above, in->reconstructedValue isn't toasted or short. */
    if (DatumGetPointer(in->reconstructedValue)) {
        reconstrValue = (kmer *) DatumGetPointer(in->reconstructedValue);
        // elog(NOTICE, "[leaf], leaf d: %s, leaf k: %d, level: %d, reconst: %s", leafValue->data, leafValue->k, level, reconstrValue->data);
    } else {
        // elog(NOTICE, "[leaf], leaf d: %s, leaf k: %d, level: %d, reconst NULL", leafValue->data, leafValue->k, level);
    }

    // Assert(reconstrValue == NULL ? level == 0 :
    //     (reconstrValue->k) == level);

    /* Reconstruct the full string represented by this leaf tuple */
    fullLen = level + (leafValue->k);
    // elog(NOTICE, "[leaf] fullLen = %d, palloc %d, level: %d", fullLen, VARHDRSZ + sizeof(int32) + fullLen + 1, level);
    if ((leafValue->k) == 0 && level > 0)
    {
        fullValue = (reconstrValue->data);
        out->leafValue = KmerPGetDatum(reconstrValue);
    }
    else
    {
        kmer       *fullKmer = (kmer *) palloc(VARHDRSZ + sizeof(int32) + fullLen + 1);
        SET_VARSIZE(fullKmer, VARHDRSZ + sizeof(int32) + fullLen + 1);
        fullKmer->k = fullLen;
        fullValue = (fullKmer->data);
        if (level)
            memcpy(fullValue, (reconstrValue->data), level);
        if ((leafValue->k) > 0)
            memcpy(fullValue + level, (char *) (leafValue->data),
                (leafValue->k));
        fullKmer->data[fullLen] = '\0';
        out->leafValue = KmerPGetDatum(fullKmer);
    }

    /* Perform the required comparison(s) */
    res = true;
    for (j = 0; j < in->nkeys; j++)
    {
        StrategyNumber strategy = in->scankeys[j].sk_strategy;
        kmer       *query = DatumGetKmerP(in->scankeys[j].sk_argument);
        qkmer      *inQkmer;
        int         queryLen = (query->k);
        int         r;
        kmer *leafKmer = DatumGetKmerP(out->leafValue);
        // elog(NOTICE, "[leaf] leafKmer: %s", leafKmer->data);

        size_t minLen = Min(queryLen, fullLen);


        switch (strategy)
        {
            case BTEqualStrategyNumber:
                r = 0;
                for (size_t i = 0; i < minLen; i++) {
                    if (!nucleotide_matches(query->data[i], fullValue[i])) {
                        r = query->data[i] - fullValue[i];
                        break;
                    }
                }
                if (r == 0)
                {
                    if (queryLen > fullLen)
                        r = -1;
                    else if (queryLen < fullLen)
                        r = 1;
                }
                res = (r == 0);
                break;
            case BTStartsWithStrategyNumber:
                r = 0;
                if (queryLen > fullLen) {
                    res = false;
                } else {
                    r = strncmp(query->data, fullValue, queryLen);
                }
                res = (r == 0);
                break;
            case BTContainsStrategyNumber:
                inQkmer = DatumGetQkmerP(in->scankeys[j].sk_argument);
                queryLen = inQkmer->k;
                minLen = Min(queryLen, fullLen);
                r = 0;
                for (size_t i = 0; i < minLen; i++) {
                    if (!nucleotide_matches(inQkmer->data[i], fullValue[i])) {
                        r = inQkmer->data[i] - fullValue[i];
                        break;
                    }
                }
                if (r == 0)
                {
                    if (queryLen > fullLen)
                        r = -1;
                    else if (queryLen < fullLen)
                        r = 1;
                }
                res = (r == 0);
                break;
            default:
                elog(ERROR, "[Leaf] unrecognized strategy number: %d",
                    in->scankeys[j].sk_strategy);
                res = false;
                break;
        }

        if (!res)
            break;              /* no need to consider remaining conditions */
    }
    PG_RETURN_BOOL(res);
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
    char *result = pstrdup(c->data);
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

Datum
qkmer_contains_swapped(PG_FUNCTION_ARGS) {
    qkmer *pattern = PG_GETARG_QKMER_P(1);
    kmer *c = PG_GETARG_KMER_P(0);
    bool result;
    if (pattern->k != c->k) {
        result = false;
    } else {
        result = contains(pattern->data, c->data);
    }
    PG_FREE_IF_COPY(pattern, 1);
    PG_FREE_IF_COPY(c, 0);   
    PG_RETURN_BOOL(result);
}