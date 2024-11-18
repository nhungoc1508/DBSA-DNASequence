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

-- **********************************s
-- * BTREE/HASH
-- **********************************
CREATE INDEX kmer_btree_idx ON kmers USING btree (kmer);
CREATE INDEX kmer_hash_idx ON kmers USING hash (kmer);