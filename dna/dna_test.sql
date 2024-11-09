create table t(id serial primary key, test dna);
insert into t Values (1,'AT');
insert into t Values (2,'ATG');
insert into t Values (3,'ATGX'); 
insert into t Values (3,'ATGC');
select test, length(test) from t;

SELECT s.id, s.sequence, k.kmer
FROM dna_sequences AS s,
     generate_kmers(s.sequence, 4) AS k(kmer);