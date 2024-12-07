#ifndef QKMER_H
#define QKMER_H

#define VALID_IUPAC_NUCLEOTIDES  16
const char iupac_nucleotides[VALID_IUPAC_NUCLEOTIDES * 2 + 1] = {'A', 'C', 'G', 'T', 'U', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N'};
/******************************************************************************
 * TYPE STRUCT
 ******************************************************************************/

typedef struct
{
    int32 k;
    char  data[MAX_KMER_LEN + 1];
} qkmer;

/******************************************************************************
 * AUXILIARY FUNCTIONS DECLARATION
 ******************************************************************************/

static qkmer *qkmer_make(int k, char *data);
static bool is_valid_qkmer(const char *str);
static bool nucleotide_matches(char iupac_code, char nucleotide);
static qkmer *qkmer_parse(const char *str);
static char *qkmer_to_str(const qkmer *c);
static bool contains(char *pattern, char *c);
static char *to_uppercase(const char *data, int length);

/******************************************************************************
 * AUXILIARY FUNCTIONS IMPLEMENTATION
 ******************************************************************************/

static qkmer *
qkmer_make(int k, char *data)
{
    qkmer *c = palloc0(sizeof(qkmer));
    c->k = k;
    strcpy(c->data, data);
    c->data[k] = '\0';
    return c;
}

static bool 
is_valid_qkmer(const char *str)
{
    for (int i = 0; str[i] != '\0'; i++) {
        bool valid = false;
        for (int j = 0; j < 16; j++) {
            if (str[i] == iupac_nucleotides[j]) {
                valid = true;
                break;
            }
        }
        if (!valid) {
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
qkmer_parse(const char *str)
{
    int k = strlen(str);
    if (k > MAX_KMER_LEN)
        ereport(ERROR, (errcode(ERRCODE_STRING_DATA_RIGHT_TRUNCATION), errmsg("Input exceeds maximum length allowed for type qkmer (32)")));
    char *upper_str = to_uppercase(str,k);
    if (!is_valid_qkmer(upper_str))
        ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION), errmsg("Invalid input syntax for type qkmer")));
    
    return qkmer_make(k, upper_str);
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
