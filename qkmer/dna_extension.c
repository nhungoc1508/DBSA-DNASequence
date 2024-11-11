#include "postgres.h"
#include "fmgr.h"
#include "utils/builtins.h"
#include "varatt.h"

PG_MODULE_MAGIC;

#define QKMER_MAX_LENGTH 32
#define VALID_IUPAC_NUCLEOTIDES "ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn"

/*****************************************************************************/

/* Structure to represent qkmer sequence */
typedef struct
{
    int32 k;
    char  *data;
} qkmer;

/* fmgr macros qkmer type */

#define DatumGetQkmerP(X)  ((qkmer *) DatumGetPointer(X))
#define QkmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_QKMER_P(n) DatumGetQkmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_QKMER_P(x) return QkmerPGetDatum(x)

/*****************************************************************************/

static qkmer *
qkmer_make(int k, char *data)
{
    qkmer *c = palloc0(sizeof(qkmer));    
    c->k = k; 
    c->data = strdup(data);
    return c;
}

static bool 
is_valid_iupac_nucleotide(char *str)
{
    for (int i=0; i<32; i++) {
        if (*str == VALID_IUPAC_NUCLEOTIDES[i]) {
            return true;
        }
    }
    return false;
}

static bool 
is_valid_qkmer(char **str)
{
    char *alt_str = *str;
    while (*alt_str) {
        if (is_valid_iupac_nucleotide(alt_str)) {
            alt_str += 1;
        } else {
            return false;
        }
    }
    return true;
}

/* Check if character matches any IUPAC nucleotide code */
static bool
nucleotide_matches(char iupac_code, char nucleotide)
{
    switch(iupac_code) // compares iupac code against cases below
    {
        case 'A': return nucleotide == 'A';
        case 'C': return nucleotide == 'C';
        case 'G': return nucleotide == 'G';
        case 'T': return nucleotide == 'T';
        case 'U': return nucleotide == 'U';
        case 'W': return nucleotide == 'A' || nucleotide == 'T';
        case 'S': return nucleotide == 'C' || nucleotide == 'G';
        case 'M': return nucleotide == 'A' || nucleotide == 'C';
        case 'K': return nucleotide == 'G' || nucleotide == 'T';
        case 'R': return nucleotide == 'A' || nucleotide == 'G';
        case 'Y': return nucleotide == 'C' || nucleotide == 'T';
        case 'B': return nucleotide == 'C' || nucleotide == 'G' || nucleotide == 'T';
        case 'D': return nucleotide == 'A' || nucleotide == 'G' || nucleotide == 'T';
        case 'H': return nucleotide == 'A' || nucleotide == 'C' || nucleotide == 'T';
        case 'V': return nucleotide == 'A' || nucleotide == 'C' || nucleotide == 'G';
        case 'N': return true;
        default: return false;
    }
}


static qkmer *
qkmer_parse(char **str)
{
    int k;
    char *data;
    bool ret = is_valid_qkmer(str);
    if (!ret)
        ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION), errmsg("Invalid input syntax for type qkmer")));
    data = strdup(*str);
    k = (int)strlen(data);
    if (k > QKMER_MAX_LENGTH)
        ereport(ERROR, (errcode(ERRCODE_STRING_DATA_RIGHT_TRUNCATION), errmsg("Input exceeds maximum length allowed for type qkmer (32)")));
    return qkmer_make(k, data);
}

static char *
qkmer_to_str(const qkmer *c) 
{
    char *result = strdup(c->data);
    return result;
}

/*****************************************************************************/

PG_FUNCTION_INFO_V1(qkmer_in);
Datum 
qkmer_in(PG_FUNCTION_ARGS) 
{
    char *str = PG_GETARG_CSTRING(0);
    PG_RETURN_QKMER_P(qkmer_parse(&str));
}

PG_FUNCTION_INFO_V1(qkmer_out);
Datum 
qkmer_out(PG_FUNCTION_ARGS) 
{
    qkmer *c = PG_GETARG_QKMER_P(0);
    char *result = strdup(c->data);
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_CSTRING(result);
}


PG_FUNCTION_INFO_V1(qkmer_cast_from_text);
Datum 
qkmer_cast_from_text(PG_FUNCTION_ARGS) 
{
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout,
                 PointerGetDatum(txt)));
    PG_RETURN_QKMER_P(qkmer_parse(&str));
}

PG_FUNCTION_INFO_V1(qkmer_cast_to_text);
Datum 
qkmer_cast_to_text(PG_FUNCTION_ARGS) 
{
    qkmer *c = PG_GETARG_QKMER_P(0);
    text *out = (text *)DirectFunctionCall1(textin,
                 PointerGetDatum(qkmer_to_str(c)));
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_TEXT_P(out);
}


/******************************************************************************
 * CONSTRUCTOR
 ******************************************************************************/
PG_FUNCTION_INFO_V1(qkmer_constructor);
Datum
qkmer_constructor(PG_FUNCTION_ARGS) 
{
    int k = PG_GETARG_INT32(0);
    char *data = text_to_cstring(PG_GETARG_TEXT_PP(1));
    PG_RETURN_QKMER_P(qkmer_make(k, data));
}

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

PG_FUNCTION_INFO_V1(qkmer_length);
Datum 
qkmer_length(PG_FUNCTION_ARGS) 
{
    qkmer *c = PG_GETARG_QKMER_P(0);
    int len = c->k;
    PG_FREE_IF_COPY(c, 0);
    PG_RETURN_INT32(len);
}
