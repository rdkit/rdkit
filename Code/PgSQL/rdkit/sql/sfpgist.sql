CREATE INDEX fpidx ON pgsfp USING gist (f);

SET rdkit.tanimoto_threshold = 0.6;
SET rdkit.dice_threshold = 0.6;


SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;

SELECT
    id, tanimoto_sml(morgan_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol, 1), f) AS sml
FROM
	pgsfp
WHERE morgan_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol, 1) % f
ORDER BY sml DESC, id;

SELECT
    id, dice_sml(morgan_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol, 1), f) AS sml
FROM
	pgsfp
WHERE morgan_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol, 1) # f
ORDER BY sml DESC, id;

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=off;

SELECT
    id, tanimoto_sml(morgan_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol, 1), f) AS sml
FROM
	pgsfp
WHERE morgan_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol, 1) % f
ORDER BY sml DESC, id;

SELECT
    id, dice_sml(morgan_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol, 1), f) AS sml
FROM
	pgsfp
WHERE morgan_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol, 1) # f
ORDER BY sml DESC, id;

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=on;

DROP INDEX fpidx;
