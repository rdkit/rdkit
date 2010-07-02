--
-- first, define the datatype.  Turn off echoing so that expected file
-- does not depend on contents of rdkit.sql.
--
SET client_min_messages = warning;
\set ECHO none
\i rdkit.sql
\set ECHO all
RESET client_min_messages;

CREATE TABLE pgmol (id int, m mol);
\copy pgmol from 'data/data'

CREATE UNIQUE INDEX mol_ididx ON pgmol (id);

SELECT count(*) FROM pgmol;

SELECT count(*) FROM pgmol WHERE m @> 'c1ccccc1';
SELECT count(*) FROM pgmol WHERE m @> 'c1cccnc1';
SELECT count(*) FROM pgmol WHERE 'c1ccccc1' <@ m;
SELECT count(*) FROM pgmol WHERE 'c1cccnc1' <@ m;

SELECT id, layered_fp(m) AS f INTO pgbfp FROM pgmol;
CREATE UNIQUE INDEX bfp_ididx ON pgbfp (id);

SELECT id, morgan_fp(m,1) AS f INTO pgsfp FROM pgmol;
CREATE UNIQUE INDEX sfp_ididx ON pgsfp (id);

SELECT id, rdkit_fp(m) AS f INTO pgavfp FROM pgmol;
SELECT id, torsion_fp(m) AS f INTO pgtorsfp FROM pgmol;
SELECT id, atompair_fp(m) AS f INTO pgpairfp FROM pgmol;

set rdkit.tanimoto_threshold=0.5;
set rdkit.dice_threshold=0.5;

SELECT 
	id, 
	tanimoto_sml(layered_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol), f) 
FROM
	 (SELECT * FROM pgbfp ORDER BY id) AS t 
WHERE layered_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol) % f  
LIMIT 10;

SELECT 
	id, 
	dice_sml(layered_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol), f) 
FROM
	 (SELECT * FROM pgbfp ORDER BY id) AS t 
WHERE layered_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol) % f  
LIMIT 10;

SELECT 
	id, 
	tanimoto_sml(layered_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol), f) 
FROM
	 (SELECT * FROM pgbfp ORDER BY id) AS t 
WHERE layered_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol) # f  
LIMIT 10;

SELECT 
	id, 
	dice_sml(layered_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol), f),
	size(f)
FROM
	 (SELECT * FROM pgbfp ORDER BY id) AS t 
WHERE layered_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol) # f  
LIMIT 10;

set rdkit.tanimoto_threshold=0.4;
SELECT 
	id, 
	tanimoto_sml(morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)))'::mol, 1), f) 
FROM
	 (SELECT * FROM pgsfp ORDER BY id) AS t 
WHERE morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)))'::mol, 1) % f  
LIMIT 10;

SELECT 
	id, 
	dice_sml(morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)))'::mol, 1), f) 
FROM
	 (SELECT * FROM pgsfp ORDER BY id) AS t 
WHERE morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)))'::mol, 1) % f  
LIMIT 10;

SELECT 
	id, 
	tanimoto_sml(morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol, 1), f) 
FROM
	 (SELECT * FROM pgsfp ORDER BY id) AS t 
WHERE morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol, 1) # f  
LIMIT 10;

SELECT 
	id, 
	dice_sml(morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol, 1), f)
FROM
	 (SELECT * FROM pgsfp ORDER BY id) AS t 
WHERE morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol, 1) # f  
LIMIT 10;


