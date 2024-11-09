CREATE EXTENSION dnasequence;

-- dna
CREATE TABLE seqs(id SERIAL PRIMARY KEY, dna dna);
INSERT INTO seqs (dna) VALUES
    ('AT'),
    ('ATG'),
    ('ATGC');
SELECT dna, length(dna) FROM seqs;

-- kmer
CREATE TABLE kmers (id SERIAL PRIMARY KEY, kmer kmer);
INSERT INTO kmers (kmer) VALUES
    ('ACGTA'),
    ('GATTC'),
    ('acgt'),
    ('gat'),
    ('AcGGtTa');
SELECT * FROM kmers;
SELECT kmer AS kmer, length(kmer) AS "length(kmer)" FROM kmers;
SELECT * FROM kmers WHERE equals('ACGTA', kmer);
SELECT * FROM kmers WHERE kmer = 'ACGTA';

-- Generate
SELECT k.kmer
FROM generate_kmers('ACGTACGT', 6) AS k(kmer);

-- SELECT s.id, s.dna, k.kmer
-- FROM seqs AS s, generate_kmers(s.dna, 4) AS k(kmer);

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