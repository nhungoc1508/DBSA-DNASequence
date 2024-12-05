CREATE EXTENSION dnasequence;
CREATE TABLE genomes(genome dna);
\copy genomes FROM './test/dna_parsed/genome_long_parsed.txt' WITH (FORMAT CSV)

CREATE TABLE sample_32mers AS
    SELECT k.kmer FROM (
        SELECT generate_kmers(genome, 32)
        FROM genomes)
    AS k(kmer);

CREATE TABLE small AS (
    SELECT * FROM sample_32mers LIMIT 1000
);

CREATE INDEX kmer_spgist_idx ON sample_32mers USING spgist(kmer);
CREATE INDEX small_idx ON small USING spgist(kmer);
SET enable_seqscan = OFF;