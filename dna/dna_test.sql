CREATE TABLE genes (
    id SERIAL PRIMARY KEY,
    sequence dna
);


INSERT INTO genes (sequence) VALUES ('ATGCGT');

SELECT * from genes;

SELECT sequence, length(sequence) from genes;