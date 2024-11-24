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

#include "data_types/dna.h"
#include "data_types/kmer.h"
#include "data_types/qkmer.h"

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

// qkmer macros
#define DatumGetQkmerP(X)  ((qkmer *) DatumGetPointer(X))
#define QkmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_QKMER_P(n) DatumGetQkmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_QKMER_P(x) return QkmerPGetDatum(x)

// SP-GiST
#define SPGIST_MAX_PREFIX_LENGTH	Max((int) (BLCKSZ - 258 * 16 - 100), 32)
#define SPG_STRATEGY_ADDITION	(10)
#define SPG_IS_COLLATION_AWARE_STRATEGY(s) ((s) > SPG_STRATEGY_ADDITION \
										 && (s) != RTPrefixStrategyNumber)

/* Struct for sorting values in picksplit */
typedef struct spgNodePtr
{
	Datum		d;
	int			i;
	int16		c;
} spgNodePtr;

#endif // DNASEQUENCE_H