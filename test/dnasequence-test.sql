CREATE EXTENSION dnasequence;

-- **********************************
-- * dna
-- **********************************
CREATE TABLE seqs(id SERIAL PRIMARY KEY, dna dna);
INSERT INTO seqs (dna) VALUES
    ('AT'),
    ('ATG'),
    ('atcg'),
    ('cgta'),
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
INSERT INTO kmers (kmer) VALUES ('A');
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
-- * IMPORT CSV DATA
-- **********************************
CREATE EXTENSION dnasequence;
CREATE TABLE genomes(genome dna);
\copy genomes FROM './test/dna_parsed/genome_long_parsed.txt' WITH (FORMAT CSV)

CREATE TABLE sample_32mers AS
    SELECT k.kmer FROM (
        SELECT generate_kmers(genome, 32)
        FROM genomes)
    AS k(kmer);

WITH kmers AS (
    SELECT generate_kmers(genome, 32) as kmer
    FROM genomes
)
SELECT kmers.kmer, count(*)
FROM kmers
GROUP BY kmers.kmer
ORDER BY count(*) DESC;

WITH kmers AS (
    SELECT generate_kmers(genome, 32) as kmer
    FROM genomes
),
kmers_count AS (
    SELECT kmers.kmer, count(*)
    FROM kmers
    GROUP BY kmers.kmer
    ORDER BY count(*) DESC
)
SELECT sum(count) AS total_count,
       count(*) AS distinct_count,
       count(*) FILTER (WHERE count = 1) AS unique_count
FROM kmers_count;

-- **********************************
-- * SP-GIST INDEX
-- **********************************
CREATE INDEX kmer_spgist_idx ON sample_32mers USING spgist(kmer);
SET enable_seqscan = OFF;
select * from sample_32mers where kmer = 'AAAGAGGCTAACAGGCTTTTGAAAAGTTATTC';
select * from sample_32mers limit 5;

-- *** Equality search ***
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE kmer = 'AAAGAGGCTAACAGGCTTTTGAAAAGTTATTC';

-- *** Prefix search ***
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE kmer ^@ 'ACG';

-- *** Pattern matching using qkmer ***
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE 'ANGTA' @> kmer;

-- **********************************
-- * SP-GIST COMPARISON
-- **********************************

-- *** Equality ***
SELECT * FROM sample_32mers WHERE kmer = 'GTGCTCCGTTCGTTCCCTCCTTGCACAGAAAG';
explain analyze select * from sample_32mers where kmer ^@ 'GTGCTCCGTTCGTTCCCTCCTTGCACA'::kmer;

-- *** Prefix search ***
drop index kmer_spgist_idx;
SET enable_seqscan = ON;
SET max_parallel_workers_per_gather = 0;
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE kmer ^@ 'AGCT';
SELECT * FROM sample_32mers WHERE kmer ^@ 'AGCT';
CREATE INDEX kmer_spgist_idx ON sample_32mers USING spgist(kmer);
SET enable_seqscan = OFF;
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE kmer ^@ 'AGCT';
SELECT * FROM sample_32mers WHERE kmer ^@ 'AGCT';

-- *** Pattern matching using qkmer ***
drop index kmer_spgist_idx;
SET enable_seqscan = ON;
SET max_parallel_workers_per_gather = 0;
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE 'ATCGAGNNNNNNNNNNNNNNNNNNNNNNNNNN' @> kmer;
SELECT * FROM sample_32mers WHERE 'ATCGAGNNNNNNNNNNNNNNNNNNNNNNNNNN' @> kmer;
CREATE INDEX kmer_spgist_idx ON sample_32mers USING spgist(kmer);
SET enable_seqscan = OFF;
EXPLAIN ANALYZE SELECT * FROM sample_32mers WHERE 'ATCGAGNNNNNNNNNNNNNNNNNNNNNNNNNN' @> kmer;
SELECT * FROM sample_32mers WHERE 'ATCGAGNNNNNNNNNNNNNNNNNNNNNNNNNN' @> kmer;

-- **********************************
-- * SYNTHETIC DATA
-- **********************************
-- * Create table where kmers will be stored
CREATE TABLE kmers_big (
    id SERIAL PRIMARY KEY,
    kmer kmer  
);

-- * Create synthetic data
DO $$
DECLARE
    total_rows BIGINT := 0;
    batch_size INT := 350000;
    max_rows BIGINT := 350000; 
BEGIN
    WHILE total_rows < max_rows LOOP
        INSERT INTO kmers_big (kmer)
        SELECT generate_kmers('ACGTACGTACGTACGTACGTACGTACGTACGT', trunc(random() * 32 + 1)::int)
        FROM generate_series(1, batch_size);

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