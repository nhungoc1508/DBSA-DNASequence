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
SELECT * FROM KMERS
WHERE equals('ACGTA', kmer)
    OR equals('acgt', kmer)
    OR equals('gAt', kmer);
-- Test = operator
SELECT * FROM kmers WHERE kmer = 'ACGTA';
SELECT * FROM kmers WHERE kmer = 'acgta';
SELECT * FROM kmers WHERE kmer = 'aCgTA';
SELECT * FROM kmers
WHERE kmer = 'ACGTA'
    OR kmer = 'acgt'
    OR kmer = 'gAt';

-- Test starts_with()
SELECT * FROM kmers WHERE starts_with('ACG', kmer);
SELECT * FROM kmers WHERE starts_with('ga', kmer);
SELECT * FROM kmers WHERE starts_with('gAT', kmer);
-- Test ^@ operator
SELECT * FROM kmers WHERE kmer ^@ 'ACG';
SELECT * FROM kmers WHERE kmer ^@ 'ga';
SELECT * FROM kmers WHERE kmer ^@ 'gAT';
SELECT * FROM kmers WHERE kmer ^@ 'tag';
SELECT * FROM kmers
WHERE kmer ^@ 'ACG'
    OR kmer ^@ 'ga'
    OR starts_with('gAT', kmer);
-- Test invalid prefix length error
SELECT * FROM kmers WHERE kmer ^@ 'ACGGTTA'; -- valid
SELECT * FROM kmers WHERE kmer ^@ 'ACGGTTAA'; -- none found