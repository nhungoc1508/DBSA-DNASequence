-- Inform users to load this file via CREATE EXTENSION
\echo Use "CREATE EXTENSION dnasequence" to load this file. \quit

/******************************************************************************
 * INPUT/OUTPUT
 ******************************************************************************/

-- ********** dna **********
CREATE OR REPLACE FUNCTION dna_in(cstring) 
    RETURNS dna 
    AS 'MODULE_PATHNAME', 'dna_in' 
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_out(dna) 
    RETURNS cstring 
    AS 'MODULE_PATHNAME', 'dna_out'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- ********** kmer **********
CREATE OR REPLACE FUNCTION kmer_in(cstring)
    RETURNS kmer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_out(kmer)
    RETURNS cstring
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/******************************************************************************
 * TYPE DEFINITIONS
 ******************************************************************************/

CREATE TYPE dna (  
    internallength = variable,
    input          = dna_in,
    output         = dna_out  
);

CREATE TYPE kmer (
    internallength = 32,
    input          = kmer_in,
    output         = kmer_out,
    -- receive        = kmer_recv,
    -- send           = kmer_send,
    alignment      = double -- ? Need to replace?
);

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/

-- ********** dna **********
CREATE FUNCTION length(dna)
    RETURNS int4
    AS 'MODULE_PATHNAME', 'dna_length'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- ********** kmer **********
CREATE OR REPLACE FUNCTION kmer(text)
    RETURNS kmer
    AS 'MODULE_PATHNAME', 'kmer_cast_from_text'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(kmer)
    RETURNS text
    AS 'MODULE_PATHNAME', 'kmer_cast_to_text'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as kmer) WITH FUNCTION kmer(text) AS IMPLICIT;
CREATE CAST (kmer as text) WITH FUNCTION text(kmer);

CREATE FUNCTION kmer(int, text)
    RETURNS kmer
    AS 'MODULE_PATHNAME', 'kmer_constructor'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION length(kmer)
    RETURNS int
    AS 'MODULE_PATHNAME', 'kmer_length'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION equals(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmer_equals'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION starts_with(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmer_starts_with'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_starts_with_swapped(kmer, kmer)
    RETURNS boolean
    AS 'SELECT starts_with($2, $1)'
    LANGUAGE SQL IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION generate_kmers(dna, integer)
    RETURNS SETOF kmer
    AS 'MODULE_PATHNAME', 'generate_kmers'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/* Additional functions for the BTree operator class */
CREATE FUNCTION kmer_lt(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmer_lt'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_le(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmer_le'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_gt(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmer_gt'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_ge(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmer_ge'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_cmp(kmer, kmer)
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmer_cmp'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/******************************************************************************
 * OPERATORS (kmer)
 ******************************************************************************/

CREATE OPERATOR = (
    LEFTARG = kmer, RIGHTARG = kmer,
    PROCEDURE = equals,
    COMMUTATOR = =
);

CREATE OPERATOR ^@ (
    LEFTARG = kmer, RIGHTARG = kmer,
    PROCEDURE = kmer_starts_with_swapped
    -- ? Commutator?
);

/* Additional operators for the BTree operator class */
CREATE OPERATOR < (
    LEFTARG = kmer, RIGHTARG = kmer,
    PROCEDURE = kmer_lt,
    COMMUTATOR = >, NEGATOR = >=
);

CREATE OPERATOR <= (
    LEFTARG = kmer, RIGHTARG = kmer,
    PROCEDURE = kmer_le,
    COMMUTATOR = >=, NEGATOR = >
);

CREATE OPERATOR > (
    LEFTARG = kmer, RIGHTARG = kmer,
    PROCEDURE = kmer_gt,
    COMMUTATOR = <, NEGATOR = <=
);

CREATE OPERATOR >= (
    LEFTARG = kmer, RIGHTARG = kmer,
    PROCEDURE = kmer_ge,
    COMMUTATOR = <=, NEGATOR = <
);

/******************************************************************************
 * OPERATOR CLASS (kmer)
 ******************************************************************************/
CREATE OPERATOR CLASS kmer_ops
    DEFAULT FOR TYPE kmer USING btree AS
        OPERATOR        1       < ,
        OPERATOR        2       <= ,
        OPERATOR        3       = ,
        OPERATOR        4       >= ,
        OPERATOR        5       > ,
        FUNCTION        1       kmer_cmp(kmer, kmer);