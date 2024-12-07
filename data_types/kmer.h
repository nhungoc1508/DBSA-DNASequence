#ifndef KMER_H
#define KMER_H

#include "common/hashfn.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>

/******************************************************************************
 * TYPE STRUCT
 ******************************************************************************/

typedef struct {
    int32 vl_len_;
    int32 k;
    char data[33];
} kmer;

/******************************************************************************
 * AUXILIARY FUNCTIONS DECLARATION
 ******************************************************************************/

static kmer *kmer_make(int k, const char *data);
static void p_whitespace(char **str);
static void ensure_end_input(char **str, bool end);
static bool is_valid_sequence(const char *sequence);
static char *to_uppercase(const char *data, int length);
static kmer *kmer_parse(const char *str);
static char *kmer_to_str(const kmer *c);
static bool starts_with(char *prefix, char *c);
static int kmer_cmp_internal(kmer *a, kmer *b);

/******************************************************************************
 * AUXILIARY FUNCTIONS IMPLEMENTATION
 ******************************************************************************/

static kmer *kmer_make(int k, const char *data) {
    kmer *c = (kmer *) palloc(VARHDRSZ + sizeof(int32) + k + 1);
    SET_VARSIZE(c, VARHDRSZ + sizeof(int32) + k + 1);
    c->k = k;
    strncpy(c->data, data, k);
    c->data[k] = '\0';
    return c;
}

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

static char *to_uppercase(const char *data, int length) {
    char *upper_str = palloc(length + 1);
    for (int i=0; i<length; i++) {
        upper_str[i] = data[i];
        if (data[i] >= 'a' && data[i] <= 'z') {
            upper_str[i] &= ~0x20;
        }
    }
    upper_str[length] = '\0';
    return upper_str;
}

static kmer *kmer_parse(const char *str) {
    int k=strlen(str);

    if (k > MAX_KMER_LEN)
        ereport(ERROR, (errcode(ERRCODE_STRING_DATA_RIGHT_TRUNCATION), errmsg("Input exceeds maximum length allowed for type kmer (32)")));
    
    char *upper_str = to_uppercase(str, k);
    if (!is_valid_sequence(upper_str))
        ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION), errmsg("Invalid input syntax for type kmer")));

    return kmer_make(k, upper_str);
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

static int kmer_cmp_internal(kmer *a, kmer *b) {
    int ret = strcmp(a->data, b->data);
    // < 0: First non-matching character in a is lower (in ASCII) than that of b
    if (ret < 0) {
        return -1;
    } else if (ret > 0) {
        return 1;
    } else {
        return 0;
    }   
}

#endif // KMER_H