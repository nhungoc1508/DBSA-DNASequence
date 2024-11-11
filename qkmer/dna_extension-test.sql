CREATE EXTENSION dna_extension

-- Create qkmer table
CREATE TABLE qkmers (
    id SERIAL PRIMARY KEY,
    qkmer qkmer
);

-- Insert valid qkmer values into table
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
INSERT INTO qkmers (qkmer) VALUES ('AATGCGTATGCTAGTACGRYSWKMACGTNNGT')

-- Test length function
SELECT length(qkmer) FROM qkmers;
