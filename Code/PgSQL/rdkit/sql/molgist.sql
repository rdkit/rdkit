CREATE INDEX molidx ON pgmol USING gist (m);

SET rdkit.tanimoto_threshold = 0.8;
SET rdkit.dice_threshold = 0.8;


SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;

SELECT count(*) FROM pgmol WHERE m @> 'c1ccccc1';
SELECT count(*) FROM pgmol WHERE m @> 'c1cccnc1';
SELECT count(*) FROM pgmol WHERE 'c1ccccc1' <@ m;
SELECT count(*) FROM pgmol WHERE 'c1cccnc1' <@ m;
SELECT count(*) FROM pgmol WHERE m @> 'c1ccccc1C(=O)N';

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=off;

SELECT count(*) FROM pgmol WHERE m @> 'c1ccccc1';
SELECT count(*) FROM pgmol WHERE m @> 'c1cccnc1';
SELECT count(*) FROM pgmol WHERE 'c1ccccc1' <@ m;
SELECT count(*) FROM pgmol WHERE 'c1cccnc1' <@ m;
SELECT count(*) FROM pgmol WHERE m @> 'c1ccccc1C(=O)N';

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=on;

DROP INDEX molidx;
