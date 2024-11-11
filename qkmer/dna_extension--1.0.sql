-- complain if script is sourced in psql, rather than via CREATE EXTENSION
\echo Use "CREATE EXTENSION dna_extension" to load this file. \quit

/******************************************************************************
 * Input/Output
 ******************************************************************************/

-- Query K-mer (qkmer) data type
CREATE OR REPLACE FUNCTION qkmer_in(cstring)
  RETURNS qkmer
  AS 'dna_extension', 'qkmer_in'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_out(qkmer)
  RETURNS cstring
  AS 'dna_extension', 'qkmer_out'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE qkmer (
  internallength = 32,
  input = qkmer_in,
  output = qkmer_out
);


CREATE OR REPLACE FUNCTION qkmer(text)
  RETURNS qkmer
  AS 'dna_extension', 'qkmer_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(qkmer)
  RETURNS text
  AS 'dna_extension', 'qkmer_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as qkmer) WITH FUNCTION qkmer(text) AS IMPLICIT;
CREATE CAST (qkmer as text) WITH FUNCTION text(qkmer);


/******************************************************************************
 * Constructor
 ******************************************************************************/

CREATE FUNCTION qkmer(int, text)
  RETURNS qkmer
  AS 'dna_extension', 'qkmer_constructor'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION length(qkmer)
  RETURNS int
  AS 'dna_extension', 'qkmer_length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;