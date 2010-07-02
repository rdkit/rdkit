CREATE INDEX fpidx ON pgfp USING gist (f);

SET rdkit.tanimoto_threshold = 0.8;
SET rdkit.dice_threshold = 0.8;


SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;

SELECT
    id, tanimoto_sml('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol::fp, f) AS sml
FROM
	pgfp
WHERE 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol::fp % f
ORDER BY sml DESC, id;

SELECT
    id, dice_sml('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol::fp, f) AS sml
FROM
	pgfp
WHERE 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol::fp # f
ORDER BY sml DESC, id;

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=off;

SELECT
    id, tanimoto_sml('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol::fp, f) AS sml
FROM
	pgfp
WHERE 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol::fp % f
ORDER BY sml DESC, id;

SELECT
    id, dice_sml('C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol::fp, f) AS sml
FROM
	pgfp
WHERE 'C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC(=C(C=C3)O)O'::mol::fp # f
ORDER BY sml DESC, id;

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=on;

DROP INDEX fpidx;
