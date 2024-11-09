CREATE EXTENSION dna;

CREATE TABLE kmers (id SERIAL, kmer kmer);
INSERT INTO kmers (kmer) VALUES
('ACGTA'),
('GATTC'),
('acgt'),
('gat'),
('AcGGtTa');

SELECT * FROM kmers;

-- Test invalid input syntax error
INSERT INTO kmers (kmer) VALUES ('ABCD');

-- Test maximum length exceeded error
INSERT INTO kmers (kmer) VALUES ('ACGTAACGTAACGTAACGTAACGTAACGTAACGTA'); -- length = 35

-- Test length()
SELECT kmer AS kmer, length(kmer) AS "length(kmer)" FROM kmers;

-- Test equals()
SELECT * FROM kmers WHERE equals('ACGTA', kmer);
SELECT * FROM kmers WHERE equals('acgta', kmer);
SELECT * FROM kmers WHERE equals('aCgTA', kmer);
SELECR * FROM
-- Test = operator
SELECT * FROM kmers WHERE kmer = 'ACGTA';
SELECT * FROM kmers WHERE kmer = 'acgta';
SELECT * FROM kmers WHERE kmer = 'aCgTA';