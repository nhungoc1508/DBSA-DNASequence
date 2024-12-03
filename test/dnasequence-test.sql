CREATE EXTENSION dnasequence;

-- **********************************
-- * dna
-- **********************************
CREATE TABLE seqs(id SERIAL PRIMARY KEY, dna dna);
INSERT INTO seqs (dna) VALUES
    ('AT'),
    ('ATG'),
    ('ATGC'),
    ('CGTA'),
    ('TAGC');

-- Test to see if values were inserted 
SELECT * FROM seqs;

-- Test validation
-- Invalid characters
INSERT INTO seqs VALUES (6,'ATGXAY'); -- Should fail

-- Test length() function
SELECT dna, length(dna) FROM seqs;

-- Generate kmers from DNA
SELECT k.kmer
FROM generate_kmers('ACGTACGT', 6) AS k(kmer);

-- **********************************
-- * kmer
-- **********************************
CREATE TABLE kmers (id SERIAL PRIMARY KEY, kmer kmer);
INSERT INTO kmers (kmer) VALUES
    ('ACGTA'),
    ('GATTC'),
    ('acgt'),
    ('gat'),
    ('AcGGtTa');

-- Test to see if values were inserted 
SELECT * FROM kmers;

-- Test invalid input syntax error
INSERT INTO kmers (kmer) VALUES ('ABCD'); -- Should fail
INSERT INTO kmers (kmer) VALUES ('12345'); -- Should fail

-- Test maximum length exceeded error
INSERT INTO kmers (kmer) VALUES ('ACGTAACGTAACGTAACGTAACGTAACGTAACGTA'); -- should fail
-- Valid edge case (length exactly 32)
INSERT INTO kmers (kmer) VALUES ('ACGTACGTACGTACGTACGTACGTACGTACGT'); -- Should pass

-- Test length()
SELECT kmer, length(kmer) FROM kmers;

-- Test equals()
SELECT * FROM kmers WHERE equals('ACGTA', kmer);
SELECT * FROM kmers WHERE kmer = 'ACGTA';
-- Test equals() with lowercase
SELECT * FROM kmers WHERE equals('acgt', kmer);
SELECT * FROM kmers WHERE kmer = 'acgt';

-- Test startswith()
SELECT * FROM kmers WHERE starts_with('ACG', kmer);
SELECT * FROM kmers WHERE kmer ^@ 'ACG';
-- Test startswith() with mixed cases
SELECT * FROM kmers WHERE starts_with('gAT', kmer);
SELECT * FROM kmers WHERE kmer ^@ 'gAT';

-- **********************************
-- * qkmer
-- **********************************
CREATE TABLE qkmers (
    id SERIAL PRIMARY KEY,
    qkmer qkmer
);

INSERT INTO qkmers (qkmer)
VALUES 
    ('AATGC'),       
    ('AAnNN'),       
    ('WgtCa'),       
    ('rACGtt'),      
    ('nTGAcGT'),     
    ('ATcGCATcG');

-- Test to see if values were inserted 
SELECT * FROM qkmers;

-- Test maximum length exceeded error (if greater than 32)
INSERT INTO qkmers (qkmer) VALUES ('AATGCGTATGCTAGTACGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN');  -- should fail

-- Test invalid nucleotide pattern matching 
-- This should raise an error since 'Z' is not a valid IUPAC code
INSERT INTO qkmers (qkmer) VALUES ('ZTGCA');  -- should fail

-- Test length of a qkmer that is exactly 32 (maximum allowed length)
INSERT INTO qkmers (qkmer) VALUES ('AATGCGTATGCTAGTACGRYSWKMACGTNNGT');

-- Test length function
SELECT qkmer, length(qkmer) FROM qkmers;

-- Test contains()
SELECT * FROM kmers WHERE contains('ANGTA', kmer);
SELECT * FROM kmers WHERE 'dhbbv'@> kmer;
SELECT * FROM kmers WHERE 'gnt'@> kmer;

-- **********************************
-- * GROUP BY
-- **********************************
SELECT k.kmer, count(*)
FROM generate_kmers('ACGTACGTGATTCACGTACGT', 5) AS k(kmer)
GROUP BY k.kmer
ORDER BY count(*) DESC;

--  kmer  | count 
----------+-------
--  TACGT |     2
--  CGTAC |     2
--  GTACG |     2
--  ACGTA |     2
--  ATTCA |     1
--  TCACG |     1
--  CACGT |     1
--  CGTGA |     1
--  ACGTG |     1
--  TTCAC |     1
--  GTGAT |     1
--  TGATT |     1
--  GATTC |     1


WITH kmers AS (
    SELECT k.kmer, count(*)
    FROM generate_kmers('ACGTACGTGATTCACGTACGT', 5) AS k(kmer)
    GROUP BY k.kmer
    ORDER BY count(*) DESC
)
SELECT sum(count) AS total_count,
       count(*) AS distinct_count,
       count(*) FILTER (WHERE count = 1) AS unique_count
FROM kmers;
--  total_count | distinct_count | unique_count 
-- -------------+----------------+--------------
--           17 |             13 |            9

-- **********************************
-- * BTREE/HASH
-- **********************************
CREATE INDEX kmer_btree_idx ON kmers USING btree (kmer);
CREATE INDEX kmer_hash_idx ON kmers USING hash (kmer);

-- **********************************
-- * IMPORT CSV DATA
-- **********************************
CREATE TABLE genomes(genome dna);
\copy genomes FROM './test/genome_parsed.txt' WITH (FORMAT CSV)

CREATE TABLE sample_32mers AS
    SELECT k.kmer FROM (
        SELECT generate_kmers(genome, 32)
        FROM genomes)
    AS k(kmer);

CREATE TABLE sample_5mers AS
    SELECT k.kmer FROM (
        SELECT generate_kmers(genome, 5)
        FROM genomes)
    AS k(kmer);

SELECT k.kmer FROM (
    SELECT generate_kmers(genome, 5)
    FROM genomes)
AS k(kmer);

SELECT k.kmer, count(*) FROM (
    SELECT generate_kmers(genome, 5)
    FROM genomes)
AS k(kmer)
GROUP BY k.kmer
ORDER BY count(*) DESC;

SELECT kmer, count(*) FROM sample_5mers GROUP BY kmer ORDER BY count(*) DESC;

-- **********************************
-- * QUERY PLANS
-- **********************************
EXPLAIN SELECT * FROM kmers; -- OK
EXPLAIN SELECT * FROM kmers WHERE kmer = 'ACGGTTA'; -- ERROR
EXPLAIN SELECT * FROM kmers WHERE kmer < 'ACGGTTA'; -- ERROR
EXPLAIN SELECT * FROM kmers WHERE equals('ACGGTTA', kmer); -- ERROR
EXPLAIN SELECT * FROM kmers WHERE kmer ^@ 'gAT'; -- ERROR
EXPLAIN SELECT * FROM kmers WHERE starts_with('ACG', kmer); -- ERROR
explain select length(kmer) from kmers; -- OK
EXPLAIN SELECT k.kmer FROM generate_kmers('ACGTACGT', 6) AS k(kmer); -- OK
EXPLAIN SELECT * FROM kmers WHERE contains('ANGTA', kmer); -- ERROR
explain select kmer(4, 'ACG'); -- OK
explain select kmer, length(kmer) from (select kmer(3, 'ACG')); -- OK
explain select kmer, text(kmer) from (select kmer(3, 'ACG')); -- OK
explain select kmer, spgist_kmer_compress(kmer) from (select kmer(3, 'ACG')); -- OK
explain select kmer, kmer_tmp_a(kmer) from kmers; -- OK
explain select kmer, kmer_tmp_b(kmer, kmer) from kmers; -- OK
explain select kmer, kmer_cmp(kmer, 'ACGGTTA') from kmers; -- OK
explain select kmer, equals(kmer, kmer) from kmers; -- OK
explain select kmer, equals(kmer, 'ACGGTTA') from kmers; -- OK
explain select * from kmers k1, kmers k2 where equals(k1.kmer, k2.kmer); -- OK
explain with k as (select 'ACGGTTA'::kmer) select * from k where kmer ^@ 'ACG'; -- OK
explain with k as (select 'ACGGTTA'::kmer) select * from k where kmer = 'ACGGTTA'; -- OK
explain select * from kmers where kmer_tmp_c(kmer, 'ACGT') is true; -- ERROR

-- **********************************
-- * SP-GIST INDEX
-- **********************************
CREATE INDEX kmer_spgist_idx ON sample_32mers USING spgist(kmer);
SET enable_seqscan = OFF;
select * from sample_32mers where kmer = 'AAAGAGGCTAACAGGCTTTTGAAAAGTTATTC';
select * from sample_5mers where kmer = 'AAAGA';
select * from sample_32mers limit 5;

-- *** Equality search ***
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE kmer = 'AAAGAGGCTAACAGGCTTTTGAAAAGTTATTC';

-- **** With SP-GIST ****
--                                                                 QUERY PLAN                                                                
-- ------------------------------------------------------------------------------------------------------------------------------------------
--  Index Only Scan using kmer_spgist_idx on sample_32mers  (cost=0.15..238.22 rows=4004 width=41) (actual time=0.050..0.051 rows=1 loops=1)
--    Index Cond: (kmer = 'AAAGAGGCTAACAGGCTTTTGAAAAGTTATTC'::kmer)
--    Heap Fetches: 0
--  Planning Time: 0.139 ms
--  Execution Time: 0.086 ms
-- (5 rows)
-- **** Without SP-GIST ****
--                                                  QUERY PLAN                                                  
-- -------------------------------------------------------------------------------------------------------------
--  Seq Scan on sample_32mers  (cost=0.00..175.11 rows=4004 width=41) (actual time=0.019..1.116 rows=1 loops=1)
--    Filter: (kmer = 'AAAGAGGCTAACAGGCTTTTGAAAAGTTATTC'::kmer)
--    Rows Removed by Filter: 8008
--  Planning Time: 0.095 ms
--  Execution Time: 1.135 ms
-- (5 rows)

-- *** Prefix search ***
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE kmer ^@ 'ACG';

-- **** With SP-GIST ****
--                                                                  QUERY PLAN                                                                 
-- --------------------------------------------------------------------------------------------------------------------------------------------
--  Index Only Scan using kmer_spgist_idx on sample_32mers  (cost=0.15..476.31 rows=2670 width=41) (actual time=2.257..2.807 rows=126 loops=1)
--    Filter: starts_with('ACG'::kmer, kmer)
--    Rows Removed by Filter: 7883
--    Heap Fetches: 0
--  Planning Time: 0.132 ms
--  Execution Time: 2.837 ms
-- (6 rows)
-- **** Without SP-GIST ****
--                                                   QUERY PLAN                                                   
-- ---------------------------------------------------------------------------------------------------------------
--  Seq Scan on sample_32mers  (cost=0.00..175.11 rows=2670 width=41) (actual time=0.088..1.271 rows=126 loops=1)
--    Filter: starts_with('ACG'::kmer, kmer)
--    Rows Removed by Filter: 7883
--  Planning Time: 0.174 ms
--  Execution Time: 1.301 ms
-- (5 rows)

-- *** Pattern matching using qkmer ***
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE 'ANGTA' @> kmer;

-- **** With SP-GIST ****
--                                                                 QUERY PLAN                                                                
-- ------------------------------------------------------------------------------------------------------------------------------------------
--  Index Only Scan using kmer_spgist_idx on sample_32mers  (cost=0.15..476.31 rows=4004 width=41) (actual time=2.149..2.150 rows=0 loops=1)
--    Filter: ('ANGTA'::qkmer @> kmer)
--    Rows Removed by Filter: 8009
--    Heap Fetches: 0
--  Planning Time: 0.090 ms
--  Execution Time: 2.171 ms
-- (6 rows)
-- **** Without SP-GIST ****
--                                                  QUERY PLAN                                                  
-- -------------------------------------------------------------------------------------------------------------
--  Seq Scan on sample_32mers  (cost=0.00..175.11 rows=4004 width=41) (actual time=0.818..0.819 rows=0 loops=1)
--    Filter: ('ANGTA'::qkmer @> kmer)
--    Rows Removed by Filter: 8009
--  Planning Time: 0.087 ms
--  Execution Time: 0.835 ms
-- (5 rows)


-- **********************************
-- * SYNTHETIC DATA
-- **********************************
-- Create table where kmers will be stored
CREATE TABLE kmers_big (
    id SERIAL PRIMARY KEY,
    kmer kmer  
);

-- Create synthetic data
/* The code inserts 1000000 random k-mers (with random lengths between 1 and 32) 
into the kmers_big table in batches of 1 million rows */

DO $$
DECLARE
    total_rows BIGINT := 0;
    batch_size INT := 1000000;
    max_rows BIGINT := 1000000; 
BEGIN
    WHILE total_rows < max_rows LOOP
        INSERT INTO kmers_big (kmer)
        SELECT generate_kmers('ACGTACGTACGTACGTACGTACGTACGTACGT', trunc(random() * 32 + 1)::int) FROM generate_series(1, batch_size);

        total_rows := total_rows + batch_size;
        RAISE NOTICE 'Inserted % rows so far', total_rows;
    END LOOP;
END;
$$;

-- Verify the inserted data
SELECT id, kmer FROM kmers_big LIMIT 10;
