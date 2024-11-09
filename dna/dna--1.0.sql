-- Inform users to load this file via CREATE EXTENSION
\echo Use "CREATE EXTENSION dna" to load this file. \quit

/****************************************************************************** 
 * Input/Output functions for the DNA type
 ******************************************************************************/ 

CREATE OR REPLACE FUNCTION dna_in(cstring) 
  RETURNS dna 
  AS 'MODULE_PATHNAME', 'dna_in' 
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_out(dna) 
  RETURNS cstring 
  AS 'MODULE_PATHNAME', 'dna_out'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- Create DNA type
CREATE TYPE dna (  
  internallength = variable,
  input          = dna_in,
  output         = dna_out  
);

CREATE FUNCTION length(dna)
  RETURNS int4
  AS 'MODULE_PATHNAME', 'dna_length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


/*------------*/
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

CREATE FUNCTION starts_with(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME', 'kmer_starts_with'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_starts_with_swapped(kmer, kmer)
    RETURNS boolean
    AS 'SELECT starts_with($2, $1)'
    LANGUAGE SQL IMMUTABLE STRICT PARALLEL SAFE;

/******************************************************************************
 * OPERATORS
 ******************************************************************************/
CREATE OPERATOR = (
    LEFTARG = kmer, RIGHTARG = kmer,
    PROCEDURE = equals,
    COMMUTATOR = =
);

CREATE OPERATOR ^@ (
    LEFTARG = kmer, RIGHTARG = kmer,
    PROCEDURE = kmer_starts_with_swapped
    -- Commutator?
);


---------

CREATE FUNCTION generate_kmers(dna, integer)
  RETURNS SETOF kmer
  AS 'MODULE_PATHNAME', 'generate_kmers'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;