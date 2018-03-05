CREATE INDEX molidx ON pgmol (m);

SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;

SELECT * FROM pgmol WHERE 
	m = 'Clc1nccnc1NN=Cc1ccc(Br)cc1';
SELECT count(*) FROM pgmol WHERE 
	m < 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O';
SELECT count(*) FROM pgmol WHERE 
	m <= 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O';
SELECT count(*) FROM pgmol WHERE 
	m = 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O';
SELECT count(*) FROM pgmol WHERE 
	m >= 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O';
SELECT count(*) FROM pgmol WHERE 
	m > 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O';

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=off;

SELECT * FROM pgmol WHERE 
	m = 'Clc1nccnc1NN=Cc1ccc(Br)cc1';
SELECT count(*) FROM pgmol WHERE 
	m < 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O';
SELECT count(*) FROM pgmol WHERE 
	m <= 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O';
SELECT count(*) FROM pgmol WHERE 
	m = 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O';
SELECT count(*) FROM pgmol WHERE 
	m >= 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O';
SELECT count(*) FROM pgmol WHERE 
	m > 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O';

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=on;

DROP INDEX molidx;

CREATE INDEX fpidx ON pgbfp (f);

SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;

SELECT * FROM pgbfp WHERE 
	f = rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);
SELECT count(*) FROM pgbfp WHERE 
	f < rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);
SELECT count(*) FROM pgbfp WHERE 
	f <= rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);
SELECT count(*) FROM pgbfp WHERE 
	f = rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);
SELECT count(*) FROM pgbfp WHERE 
	f >= rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);
SELECT count(*) FROM pgbfp WHERE 
	f > rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=off;

SELECT * FROM pgbfp WHERE 
	f = rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);
SELECT count(*) FROM pgbfp WHERE 
	f < rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);
SELECT count(*) FROM pgbfp WHERE 
	f <= rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);
SELECT count(*) FROM pgbfp WHERE 
	f = rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);
SELECT count(*) FROM pgbfp WHERE 
	f >= rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);
SELECT count(*) FROM pgbfp WHERE 
	f > rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol);

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=on;

DROP INDEX fpidx;

CREATE INDEX fpidx ON pgsfp (f);

SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;

SELECT * FROM pgsfp WHERE 
	f = morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);
SELECT count(*) FROM pgsfp WHERE 
	f < morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);
SELECT count(*) FROM pgsfp WHERE 
	f <= morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);
SELECT count(*) FROM pgsfp WHERE 
	f = morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);
SELECT count(*) FROM pgsfp WHERE 
	f >= morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);
SELECT count(*) FROM pgsfp WHERE 
	f > morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=off;

SELECT * FROM pgsfp WHERE 
	f = morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);
SELECT count(*) FROM pgsfp WHERE 
	f < morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);
SELECT count(*) FROM pgsfp WHERE 
	f <= morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);
SELECT count(*) FROM pgsfp WHERE 
	f = morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);
SELECT count(*) FROM pgsfp WHERE 
	f >= morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);
SELECT count(*) FROM pgsfp WHERE 
	f > morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol, 1);

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=on;

DROP INDEX fpidx;
