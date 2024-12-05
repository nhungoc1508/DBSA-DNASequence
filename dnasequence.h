#ifndef DNASEQUENCE_H
#define DNASEQUENCE_H

#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"

#include "access/spgist.h"
#include "common/int.h"
#include "utils/datum.h"
#include "utils/pg_locale.h"
#include "utils/varlena.h"

#include "c.h"

#include "data_types/dna.h"
#include "data_types/kmer.h"
#include "data_types/qkmer.h"

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

// qkmer macros
#define DatumGetQkmerP(X)  ((qkmer *) DatumGetPointer(X))
#define QkmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_QKMER_P(n) DatumGetQkmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_QKMER_P(x) return QkmerPGetDatum(x)

// SP-GiST

/* Struct for sorting values in picksplit */
typedef struct spgNodePtr
{
	Datum		d;
	int			i;
	int16		c;
} spgNodePtr;

typedef uint16 StrategyNumber;
#define BTEqualStrategyNumber           1
#define BTStartsWithStrategyNumber		2
#define BTContainsStrategyNumber		3

#endif // DNASEQUENCE_H