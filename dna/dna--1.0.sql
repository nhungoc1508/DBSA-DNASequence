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