#ifndef DNASEQUENCE_H
#define DNASEQUENCE_H

#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"

#include "data_types/dna.h"
#include "data_types/kmer.h"

// PG_MODULE_MAGIC;

// dna macros
#define DatumGetDnaP(X) ((dna *) PG_DETOAST_DATUM(X))
#define DnaPGetDatum(X) PointerGetDatum(X)
#define PG_GETARG_DNA_P(n) DatumGetDnaP(PG_GETARG_DATUM(n))
#define PG_RETURN_DNA_P(x) return DnaPGetDatum(x)

// kmer macros
#define DatumGetKmerP(X) ((kmer *) DatumGetPointer(X))
#define KmerPGetDatum(X) PointerGetDatum(X)
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x)

#endif // DNASEQUENCE_H