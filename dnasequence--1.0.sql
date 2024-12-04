-- Inform users to load this file via CREATE EXTENSION
\echo Use "CREATE EXTENSION dnasequence" to load this file. \quit

/******************************************************************************
 * INPUT/OUTPUT
 ******************************************************************************/

-- ********** dna **********
CREATE OR REPLACE FUNCTION dna_in(cstring) 
    RETURNS dna 
    AS 'MODULE_PATHNAME', 'dna_in' 
    LANGUAGE C IMMUTABLE;

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

-- ********** qkmer **********
CREATE OR REPLACE FUNCTION qkmer_in(cstring)
    RETURNS qkmer
    AS 'MODULE_PATHNAME', 'qkmer_in'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_out(qkmer)
    RETURNS cstring
    AS 'MODULE_PATHNAME', 'qkmer_out'
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
    internallength = variable,
    input          = kmer_in,
    output         = kmer_out,
    alignment      = double
);

CREATE TYPE qkmer (
    internallength = 32,
    input = qkmer_in,
    output = qkmer_out
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

/* Additional functions for the hash operator class */
CREATE FUNCTION kmer_hash(kmer)
    RETURNS integer
    AS 'MODULE_PATHNAME', 'kmer_hash'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/* Additional functions for the SP-GiST operator class */
CREATE FUNCTION spgist_kmer_config(internal, internal)
    RETURNS void
    AS 'MODULE_PATHNAME', 'spgist_kmer_config'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION spgist_kmer_choose(internal, internal)
    RETURNS void
    AS 'MODULE_PATHNAME', 'spgist_kmer_choose'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION spgist_kmer_picksplit(internal, internal)
    RETURNS void
    AS 'MODULE_PATHNAME', 'spgist_kmer_picksplit'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION spgist_kmer_inner_consistent(internal, internal)
    RETURNS void
    AS 'MODULE_PATHNAME', 'spgist_kmer_inner_consistent'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION spgist_kmer_leaf_consistent(internal, internal)
    RETURNS bool
    AS 'MODULE_PATHNAME', 'spgist_kmer_leaf_consistent'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- CREATE FUNCTION spgist_kmer_compress(kmer)
--     RETURNS text
--     AS 'MODULE_PATHNAME', 'spgist_kmer_compress'
--     LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- ********** qkmer **********
CREATE OR REPLACE FUNCTION qkmer(text)
    RETURNS qkmer
    AS 'MODULE_PATHNAME', 'qkmer_cast_from_text'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(qkmer)
    RETURNS text
    AS 'MODULE_PATHNAME', 'qkmer_cast_to_text'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as qkmer) WITH FUNCTION qkmer(text) AS IMPLICIT;
CREATE CAST (qkmer as text) WITH FUNCTION text(qkmer);

CREATE FUNCTION qkmer(int, text)
    RETURNS qkmer
    AS 'MODULE_PATHNAME', 'qkmer_constructor'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION length(qkmer)
    RETURNS int
    AS 'MODULE_PATHNAME', 'qkmer_length'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION contains(qkmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'qkmer_contains'
    LANGUAGE C IMMUTABLE STRICT PARALLEL RESTRICTED;

/******************************************************************************
 * OPERATORS (kmer)
 ******************************************************************************/

CREATE OPERATOR = (
    LEFTARG = kmer, RIGHTARG = kmer,
    PROCEDURE = equals,
    COMMUTATOR = =,
    HASHES
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
 * OPERATORS (qkmer)
 ******************************************************************************/

CREATE OPERATOR @> (
    LEFTARG = qkmer, RIGHTARG = kmer,
    PROCEDURE = contains
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

CREATE OPERATOR CLASS kmer_hash_ops
    DEFAULT FOR TYPE kmer USING hash AS
        OPERATOR        1       = (kmer, kmer),
        FUNCTION        1       kmer_hash(kmer);

CREATE OPERATOR CLASS kmer_spgist_ops
    DEFAULT FOR TYPE kmer USING spgist AS
        OPERATOR        1       < ,
        OPERATOR        2       <= ,
        OPERATOR        3       = ,
        OPERATOR        4       >= ,
        OPERATOR        5       > ,
        FUNCTION        1       spgist_kmer_config(internal, internal),
        FUNCTION        2       spgist_kmer_choose(internal, internal),
        FUNCTION        3       spgist_kmer_picksplit(internal, internal),
        FUNCTION        4       spgist_kmer_inner_consistent(internal, internal),
        FUNCTION        5       spgist_kmer_leaf_consistent(internal, internal),
        -- FUNCTION        6       spgist_kmer_compress(kmer),
        STORAGE         kmer;
