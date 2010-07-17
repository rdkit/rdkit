CREATE INDEX fpidx ON pgbfp USING gist (f);

SET rdkit.tanimoto_threshold = 0.7;
SET rdkit.dice_threshold = 0.8;


SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;

SELECT
    id, tanimoto_sml(layered_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol), f) AS sml
FROM
	pgbfp
WHERE layered_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol) % f
ORDER BY sml DESC, id limit 10;

SELECT
    id, dice_sml(layered_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol), f) AS sml
FROM
	pgbfp
WHERE layered_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol) # f
ORDER BY sml DESC, id;

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=off;

SELECT
    id, tanimoto_sml(layered_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol), f) AS sml
FROM
	pgbfp
WHERE layered_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol) % f
ORDER BY sml DESC, id limit 10;

SELECT
    id, dice_sml(layered_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol), f) AS sml
FROM
	pgbfp
WHERE layered_fp('O=C1CC(OC2=CC=CC=C12)C1=CC=CC=C1'::mol) # f
ORDER BY sml DESC, id limit 10;

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=on;

DROP INDEX fpidx;
