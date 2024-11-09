-- complain if script is sourced in psql, rather than via CREATE EXTENSION
\echo Use "CREATE EXTENSION dna" to load this file. \quit

/******************************************************************************
 * Input/Output
 ******************************************************************************/

CREATE OR REPLACE FUNCTION kmer_in(cstring)
    RETURNS kmer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_out(kmer)
    RETURNS cstring
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- CREATE OR REPLACE FUNCTION kmer_recv(internal)
--     RETURNS kmer
--     AS 'MODULE_PATHNAME'
--     LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- CREATE OR REPLACE FUNCTION kmer_send(kmer)
--     RETURNS bytea
--     AS 'MODULE_PATHNAME'
--     LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE kmer (
    internallength = 32,
    input          = kmer_in,
    output         = kmer_out,
    -- receive        = kmer_recv,
    -- send           = kmer_send,
    alignment      = double -- ? Need to replace?
);

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

/******************************************************************************
 * CONSTRUCTOR
 ******************************************************************************/

CREATE FUNCTION kmer(int, text)
    RETURNS kmer
    AS 'MODULE_PATHNAME', 'kmer_constructor'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/******************************************************************************
 * FUNCTIONS
 ******************************************************************************/
CREATE FUNCTION length(kmer)
    RETURNS int
    AS 'MODULE_PATHNAME', 'kmer_length'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION equals(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmer_equals'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/******************************************************************************
 * OPERATORS
 ******************************************************************************/
CREATE OPERATOR = (
    LEFTARG = kmer, RIGHTARG = kmer,
    PROCEDURE = equals,
    COMMUTATOR = =
);