CREATE EXTENSION dnasequence;

-- **********************************
-- * dna
-- **********************************
CREATE TABLE seqs(id SERIAL PRIMARY KEY, dna dna);
INSERT INTO seqs (dna) VALUES
    ('AT'),
    ('ATG'),
    ('ATGC');
-- Test validation
INSERT INTO seqs VALUES (4,'ATGX');
SELECT dna, length(dna) FROM seqs;
-- Generate
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
-- Test invalid input syntax error
INSERT INTO kmers (kmer) VALUES ('ABCD');
-- Test maximum length exceeded error
INSERT INTO kmers (kmer) VALUES ('ACGTAACGTAACGTAACGTAACGTAACGTAACGTA'); -- length = 35
-- Test length()
SELECT kmer, length(kmer) FROM kmers;
-- Test equals()
SELECT * FROM kmers WHERE equals('ACGTA', kmer);
SELECT * FROM kmers WHERE kmer = 'ACGTA';
-- Test startswith()
SELECT * FROM kmers WHERE starts_with('ACG', kmer);
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
    ('AANNN'),       
    ('WGTCA'),       
    ('RACGTT'),      
    ('NTGACGT'),     
    ('ATGGCATCG');
-- Test to see if values were inserted 
SELECT * FROM qkmers;
-- Test maximum length exceeded error (if greater than 32)
INSERT INTO qkmers (qkmer) VALUES ('AATGCGTATGCTAGTACGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN');  
-- Test invalid nucleotide pattern matching 
-- This should raise an error since 'Z' is not a valid IUPAC code
INSERT INTO qkmers (qkmer) VALUES ('ZTGCA');  
-- Test length of a qkmer that is exactly 32 (maximum allowed length), should work
INSERT INTO qkmers (qkmer) VALUES ('AATGCGTATGCTAGTACGRYSWKMACGTNNGT');
-- Test length function
SELECT qkmer, length(qkmer) FROM qkmers;
-- Test contains()
SELECT * FROM kmers WHERE contains('ANGTA', kmer);
SELECT * FROM kmers WHERE 'dhbbv'@> kmer;

-- **********************************
-- * GROUP BY
-- **********************************
SELECT k.kmer, count(*)
FROM generate_kmers('ACGTACGTGATTCACGTACGT', 5) AS k(kmer)
GROUP BY k.kmer
ORDER BY count(*) DESC;

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
CREATE INDEX kmer_spgist_idx ON sample_5mers USING spgist(kmer);
SET enable_seqscan = OFF;

-- *** Equality search ***
EXPLAIN ANALYZE SELECT * FROM sample_5mers WHERE kmer = 'AAGAA';
-- **** With SP-GIST ****
--                                                                QUERY PLAN                                                               
-- ----------------------------------------------------------------------------------------------------------------------------------------
--  Index Only Scan using kmer_spgist_idx on sample_5mers  (cost=0.15..83.36 rows=2012 width=14) (actual time=0.566..0.570 rows=19 loops=1)
--    Index Cond: (kmer = 'AAGAA'::kmer)
--    Heap Fetches: 0
--  Planning Time: 0.129 ms
--  Execution Time: 0.955 ms
-- **** Without SP-GIST ****
--                                                           QUERY PLAN                                                          
-- ------------------------------------------------------------------------------------------------------------------------------
--  Seq Scan on sample_5mers  (cost=10000000000.00..10000000072.29 rows=2012 width=14) (actual time=0.032..0.575 rows=21 loops=1)
--    Filter: (kmer = 'AAGAA'::kmer)
--    Rows Removed by Filter: 4002
--  Planning Time: 0.106 ms
--  Execution Time: 0.603 ms
-- (5 rows)

-- *** Prefix search ***
EXPLAIN ANALYZE SELECT * FROM sample_5mers WHERE kmer ^@ 'ACG';
-- **** With SP-GIST ****
--                                                                QUERY PLAN                                                                
-- -----------------------------------------------------------------------------------------------------------------------------------------
--  Index Only Scan using kmer_spgist_idx on sample_5mers  (cost=0.15..158.55 rows=1341 width=14) (actual time=1.323..1.698 rows=28 loops=1)
--    Filter: starts_with('ACG'::kmer, kmer)
--    Rows Removed by Filter: 3995
--    Heap Fetches: 0
--  Planning Time: 0.179 ms
--  Execution Time: 1.752 ms
-- (6 rows)
-- **** Without SP-GIST ****
--                                                           QUERY PLAN                                                          
-- ------------------------------------------------------------------------------------------------------------------------------
--  Seq Scan on sample_5mers  (cost=10000000000.00..10000000072.29 rows=1341 width=14) (actual time=0.069..0.599 rows=30 loops=1)
--    Filter: starts_with('ACG'::kmer, kmer)
--    Rows Removed by Filter: 3993
--  Planning Time: 0.156 ms
--  Execution Time: 0.626 ms

-- *** Pattern matching using qkmer ***
EXPLAIN ANALYZE SELECT * FROM sample_5mers WHERE 'ANGTA' @> kmer;
-- **** With SP-GIST ****
--                                                                QUERY PLAN                                                               
-- ----------------------------------------------------------------------------------------------------------------------------------------
--  Index Only Scan using kmer_spgist_idx on sample_5mers  (cost=0.15..158.55 rows=2012 width=14) (actual time=0.656..1.585 rows=5 loops=1)
--    Filter: ('ANGTA'::qkmer @> kmer)
--    Rows Removed by Filter: 4018
--    Heap Fetches: 0
--  Planning Time: 0.182 ms
--  Execution Time: 1.637 ms
-- (6 rows)
-- **** Without SP-GIST ****
--                                                          QUERY PLAN                                                          
-- -----------------------------------------------------------------------------------------------------------------------------
--  Seq Scan on sample_5mers  (cost=10000000000.00..10000000072.29 rows=2012 width=14) (actual time=0.187..0.680 rows=5 loops=1)
--    Filter: ('ANGTA'::qkmer @> kmer)
--    Rows Removed by Filter: 4018
--  Planning Time: 0.125 ms
--  Execution Time: 0.706 ms
-- (5 rows)