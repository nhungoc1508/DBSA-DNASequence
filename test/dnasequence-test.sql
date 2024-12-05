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
-- Insert a valid DNA sequence with only one character
INSERT INTO seqs (dna) VALUES ('A'); 
-- Insert a DNA sequence with repeated valid characters
INSERT INTO seqs (dna) VALUES ('AAAAAAA'); 
-- Insert very long sequences
INSERT INTO seqs (dna) VALUES (('ACGT' || repeat('GT', 3500))::text::dna);


-- Invalid characters
INSERT INTO seqs VALUES (6,'ATGXAY'); -- Should fail
INSERT INTO seqs (dna) VALUES (''); -- Should fail 

-- Test length() function
SELECT dna, length(dna) FROM seqs;
SELECT dna, length(dna) FROM seqs WHERE length(dna) = 1;
SELECT length(dna) FROM seqs WHERE length(dna) > 5000;



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

-- Test validation 
-- Insert only one character
INSERT INTO kmers (kmer) VALUES ('A')
-- Valid edge case (length exactly 32)
INSERT INTO kmers (kmer) VALUES ('ACGTACGTACGTACGTACGTACGTACGTACGT'); -- Should pass

-- Test invalid input syntax error
INSERT INTO kmers (kmer) VALUES ('ABCD'); -- Should fail
INSERT INTO kmers (kmer) VALUES ('12345'); -- Should fail
INSERT INTO kmers (kmer) VALUES ('ACGTX'); -- Should fail due to 'X'
INSERT INTO kmers (kmer) VALUES (''); -- Should fail 
-- Test maximum length exceeded error
INSERT INTO kmers (kmer) VALUES ('ACGTAACGTAACGTAACGTAACGTAACGTAACGTA'); -- should fail


-- Test length()
SELECT kmer, length(kmer) FROM kmers;
SELECT kmer, length(kmer) FROM kmers WHERE length(kmer) BETWEEN 1 AND 32;


-- Test equals()
SELECT * FROM kmers WHERE equals('ACGTA', kmer);
SELECT * FROM kmers WHERE kmer = 'ACGTA';
-- Test equals() with lowercase
SELECT * FROM kmers WHERE equals('acgt', kmer);
SELECT * FROM kmers WHERE kmer = 'acgt';
-- Test invalid codes
SELECT * FROM kmers WHERE equals('ACGTX', kmer); -- Should fail 
SELECT * FROM kmers WHERE kmer = '12345'; -- Should fail 


-- Test startswith()
SELECT * FROM kmers WHERE starts_with('ACG', kmer);
SELECT * FROM kmers WHERE kmer ^@ 'ACG';
-- Test startswith() with mixed cases
SELECT * FROM kmers WHERE starts_with('gAT', kmer);
SELECT * FROM kmers WHERE kmer ^@ 'gAT';
-- Test invalid codes
SELECT * FROM kmers WHERE starts_with('ACGTX', kmer); -- Should fail 
SELECT * FROM kmers WHERE kmer ^@ '123'; -- Should fail 


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

-- Test validation
-- Test length of a qkmer that is exactly 32 (maximum allowed length)
INSERT INTO qkmers (qkmer) VALUES ('AATGCGTATGCTAGTACGRYSWKMACGTNNGT');


-- Test invalid nucleotides
INSERT INTO qkmers (qkmer) VALUES ('ZTGCA');  -- should fail
-- Test maximum length exceeded error (if greater than 32)
INSERT INTO qkmers (qkmer) VALUES ('AATGCGTATGCTAGTACGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN');  -- should fail
INSERT INTO qkmers (qkmer) VALUES (''); -- Should fail 

-- Test length function
SELECT qkmer, length(qkmer) FROM qkmers;
SELECT kmer, length(kmer) FROM kmers WHERE length(kmer) BETWEEN 1 AND 32;


-- Test contains()
SELECT * FROM kmers WHERE contains('ANGTA', kmer);
SELECT * FROM kmers WHERE 'dhbbv'@> kmer;
SELECT * FROM kmers WHERE 'gnt'@> kmer;
-- Test invalid codes
SELECT * FROM kmers WHERE contains('ACGTX', kmer); -- Should fail
SELECT * FROM kmers WHERE 'QWDP'@> kmer;

-- Test contains and equals to
SELECT * FROM kmers WHERE contains('ACGTA', kmer) AND equals('ACGTA', kmer);
-- Test contains and starts with
SELECT * FROM kmers WHERE contains('NNGT', kmer) AND starts_with('ACG', kmer); 

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


-- Test with invalid characters -- should fail
SELECT k.kmer, count(*)
FROM generate_kmers('AOISJCGTGATTSADFDCACGTACAFGT', 5) AS k(kmer)
GROUP BY k.kmer
ORDER BY count(*) DESC;
-- ERROR:  Invalid character in DNA sequence
-- LINE 2: FROM generate_kmers('AOISJCGTGATTSADFDCACGTACAFGT', 5) AS k(...
--                             ^

-- Test with invalid length (33) -- should fail
SELECT k.kmer, count(*)
FROM generate_kmers('ACGTACGTGATTCACGTACGT', 33) AS k(kmer)
GROUP BY k.kmer
ORDER BY count(*) DESC;

-- Test with invalid length (0) -- should fail
SELECT k.kmer, count(*)
FROM generate_kmers('ACGTACGTGATTCACGTACGT', 0) AS k(kmer)
GROUP BY k.kmer
ORDER BY count(*) DESC;
-- ERROR:  Invalid k value: must be between 1 and min(sequence length, 32)

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
\copy genomes FROM './test/dna_parsed/genome_long_parsed.txt' WITH (FORMAT CSV)

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
-- * SP-GIST COMPARISON
-- **********************************
CREATE INDEX kmer_spgist_idx ON sample_32mers USING spgist(kmer);

-- *** Equality ***
SELECT * FROM sample_32mers WHERE kmer = 'GTGCTCCGTTCGTTCCCTCCTTGCACAGAAAG';
SELECT * FROM small WHERE kmer = 'GTGCTCCGTTCGTTCCCTCCTTGCACAGAAAG';
--                                                     QUERY PLAN                                                    
-- ------------------------------------------------------------------------------------------------------------------
--  Bitmap Heap Scan on small  (cost=32.02..54.27 rows=500 width=41) (actual time=0.643..0.646 rows=1 loops=1)
--    Recheck Cond: (kmer = 'GTGCTCCGTTCGTTCCCTCCTTGCACAGAAAG'::kmer)
--    Heap Blocks: exact=1
--    ->  Bitmap Index Scan on sm_idx  (cost=0.00..31.89 rows=500 width=0) (actual time=0.608..0.609 rows=1 loops=1)
--          Index Cond: (kmer = 'GTGCTCCGTTCGTTCCCTCCTTGCACAGAAAG'::kmer)
--  Planning Time: 0.244 ms
--  Execution Time: 0.745 ms
-- (7 rows)

explain analyze select * from small where kmer ^@ 'GTGCTCCGTTCGTTCCCTCCTTGCACA'::kmer;
--                                                       QUERY PLAN                                                      
-- ----------------------------------------------------------------------------------------------------------------------
--  Bitmap Heap Scan on small  (cost=61.23..89.73 rows=333 width=41) (actual time=7.420..8.096 rows=1 loops=1)
--    Filter: starts_with('GTGCTCCGTTCGTTCCCTCCTTGCACA'::kmer, kmer)
--    Rows Removed by Filter: 999
--    Heap Blocks: exact=10
--    ->  Bitmap Index Scan on sm_idx  (cost=0.00..61.14 rows=1000 width=0) (actual time=7.284..7.285 rows=1000 loops=1)
--  Planning Time: 0.425 ms
--  Execution Time: 8.218 ms
-- (7 rows)

-- *** Prefix search ***
drop index kmer_spgist_idx;
SET enable_seqscan = ON;
SET max_parallel_workers_per_gather = 0;
\o test/prefix_search/seqscan.txt
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE kmer ^@ 'AGCT';
SELECT * FROM sample_32mers WHERE kmer ^@ 'AGCT';
\o
CREATE INDEX kmer_spgist_idx ON sample_32mers USING spgist(kmer);
SET enable_seqscan = OFF;
\o test/prefix_search/indexscan.txt
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE kmer ^@ 'AGCT';
SELECT * FROM sample_32mers WHERE kmer ^@ 'AGCT';
\o

select * from sample_32mers where kmer ^@ 'AGCTAAGAAACAACTCGCTCTGCACAGG';

-- *** Pattern matching using qkmer ***
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE 'AAAGAGNNNNNNNNNNNNNNNNNNNNNNNNNN' @> kmer;

drop index kmer_spgist_idx;
SET enable_seqscan = ON;
SET max_parallel_workers_per_gather = 0;
\o test/pattern_matching/seqscan.txt
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE 'ATCGAGNNNNNNNNNNNNNNNNNNNNNNNNNN' @> kmer;
SELECT * FROM sample_32mers WHERE 'ATCGAGNNNNNNNNNNNNNNNNNNNNNNNNNN' @> kmer;
\o
CREATE INDEX kmer_spgist_idx ON sample_32mers USING spgist(kmer);
SET enable_seqscan = OFF;
\o test/pattern_matching/indexscan.txt
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE 'ATCGAGNNNNNNNNNNNNNNNNNNNNNNNNNN' @> kmer;
SELECT * FROM sample_32mers WHERE 'ATCGAGNNNNNNNNNNNNNNNNNNNNNNNNNN' @> kmer;
\o

-- * SYNTHETIC DATA
-- **********************************
-- Create table where kmers will be stored
CREATE TABLE kmers_big (
    id SERIAL PRIMARY KEY,
    kmer kmer  
);

-- Create synthetic data
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
SELECT COUNT(*) FROM kmers_big;
SELECT id, kmer FROM kmers_big LIMIT 10;

--- Test kmers_big
SELECT id, kmer, length(kmer) FROM kmers_big LIMIT 20;
SELECT count(*) FROM kmers_big WHERE length(kmer) > 32; -- ensure no kmer is greater than 32
SELECT id, kmer, length(kmer) FROM kmers_big WHERE length(kmer) = 32 LIMIT 20;
SELECT count(*) FROM kmers_big WHERE equals('ACGTACGT', kmer);
SELECT count(*) FROM kmers_big WHERE starts_with('ACG', kmer);
SELECT count(*) FROM kmers_big WHERE contains('NGT', kmer);

-- Count distinct kmers 
SELECT kmer, COUNT(*) 
FROM kmers_big
GROUP BY kmer
ORDER BY COUNT(*) DESC
LIMIT 10;
-- Check distribution of kmers_big
SELECT length(kmer) AS kmer_length, COUNT(*) 
FROM kmers_big
GROUP BY length(kmer)
ORDER BY length(kmer);
-- Check unique kmers
SELECT DISTINCT kmer FROM kmers_big ORDER BY kmer LIMIT 10;
-- Check if there are empty kmers
SELECT * FROM kmers_big WHERE kmer = '';


-- Performance Test
EXPLAIN ANALYZE SELECT * FROM kmers_big WHERE kmer = 'ACGTTGCA';
EXPLAIN ANALYZE SELECT * FROM kmers_big WHERE kmer ^@ 'ACG';


