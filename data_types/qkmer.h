#ifndef QKMER_H
#define QKMER_H

#define QKMER_MAX_LENGTH 32
#define VALID_IUPAC_NUCLEOTIDES "ACGTUWSMKRYBDHVNacgtuwsmkrybdhvn"

/******************************************************************************
 * TYPE STRUCT
 ******************************************************************************/

typedef struct
{
    int32 k;
    char  data[32];
} qkmer;

/******************************************************************************
 * AUXILIARY FUNCTIONS DECLARATION
 ******************************************************************************/

static qkmer *qkmer_make(int k, char *data);
static bool is_valid_iupac_nucleotide(char *str);
static bool is_valid_qkmer(char **str);
static bool nucleotide_matches(char iupac_code, char nucleotide);
static qkmer *qkmer_parse(char **str);
static char *qkmer_to_str(const qkmer *c);
static bool contains(char *pattern, char *c);

/******************************************************************************
 * AUXILIARY FUNCTIONS IMPLEMENTATION
 ******************************************************************************/

static qkmer *
qkmer_make(int k, char *data)
{
    // qkmer *c = palloc0(sizeof(qkmer)); 
    qkmer *c = (qkmer *) palloc(VARHDRSZ + k + 1);
    SET_VARSIZE(c, VARHDRSZ + k + 1);
    c->k = k;
    strcpy(c->data, data);
    c->data[k] = '\0';
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
    iupac_code &= ~0x20;
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

static bool contains(char *pattern, char *c) {

    int len = strlen(pattern);
    for (int i = 0; i < len; i++) {
        if (!nucleotide_matches(pattern[i], c[i])) {
            return false;
        }
    }
    return true;
}

#endif // QKMER_H
