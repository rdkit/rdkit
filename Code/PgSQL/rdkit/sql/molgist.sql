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

-- ###############################
-- github issue #525
CREATE TABLE chemblmol (molregno int, m mol);
\copy chemblmol from 'data/chembl20_100.csv' (format csv)
CREATE INDEX mol_idx2 ON chemblmol using gist(m);
select * into chemblmol2 from chemblmol order by molregno asc limit 10;
CREATE INDEX mol_idx22 ON chemblmol2 using gist(m);


-- start with a direct seq scan to verify that there is a result
SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;
select count(*) from chemblmol where m@='Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl'::mol;
select count(*) from chemblmol join chemblmol2 using (m);
set rdkit.do_chiral_sss=true;
select count(*) from chemblmol c2 join chemblmol using (m);
set rdkit.do_chiral_sss=false;
select count(*) from chemblmol c2 join chemblmol using (m);

-- now enable the index to trigger the bug:
SET enable_indexscan=on;
SET enable_bitmapscan=off;
SET enable_seqscan=on;
select count(*) from chemblmol where m@='Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl'::mol;
select count(*) from chemblmol join chemblmol2 using (m);
set rdkit.do_chiral_sss=true;
select count(*) from chemblmol c2 join chemblmol using (m);
set rdkit.do_chiral_sss=false;
select count(*) from chemblmol c2 join chemblmol using (m);


-- build index on qmols
CREATE TABLE pgqmol (id int, m qmol);
\copy pgqmol from 'data/qmol_data'
CREATE INDEX molidx ON pgqmol USING gist (m);

SET rdkit.tanimoto_threshold = 0.8;
SET rdkit.dice_threshold = 0.8;


SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;

SELECT count(*) FROM pgqmol WHERE m <@ 'c1ccccc1';
SELECT count(*) FROM pgqmol WHERE m <@ 'c1cccnc1';
SELECT count(*) FROM pgqmol WHERE 'c1ccccc1' @> m;
SELECT count(*) FROM pgqmol WHERE 'c1cccnc1' @> m;
SELECT count(*) FROM pgqmol WHERE m<@ 'c1ccccc1C(=O)N';

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=off;

SELECT count(*) FROM pgqmol WHERE m <@ 'c1ccccc1';
SELECT count(*) FROM pgqmol WHERE m <@ 'c1cccnc1';
SELECT count(*) FROM pgqmol WHERE 'c1ccccc1' @> m;
SELECT count(*) FROM pgqmol WHERE 'c1cccnc1' @> m;
SELECT count(*) FROM pgqmol WHERE m<@ 'c1ccccc1C(=O)N';

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=on;

DROP INDEX molidx;
