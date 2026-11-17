-- Nearest-neighbour search (ORDER BY <%> / <#>) through a bfp GiST index
-- must return the same top 10 as a plain scan.
-- The index needs more than one leaf page for this to be a real test, so we
-- use 200 similar benzamides with 512-bit Morgan fingerprints.
SET extra_float_digits=0;
SET rdkit.morgan_fp_size = 512;

CREATE TABLE knn_mols (id int, smiles text);
\copy knn_mols from 'data/bfpknn_data'

CREATE TABLE knn_fps AS
SELECT id, morganbv_fp(mol_from_smiles(smiles::cstring)) AS fp FROM knn_mols;
CREATE INDEX knn_fps_idx ON knn_fps USING gist (fp);

-- at least 3 pages = a root page plus two or more leaf pages
SELECT pg_relation_size('knn_fps_idx') / current_setting('block_size')::int >= 3 AS multi_page;

-- the query molecule is not in the table
CREATE TABLE knn_query AS
SELECT morganbv_fp('CCc1ccc(C(=O)NCc2ccccc2)cc1'::mol) AS fp;


-- Tanimoto

-- exact answer: compute every distance and sort, no index
SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;
CREATE TABLE knn_exact AS
SELECT id, 1.0 - tanimoto_sml(fp, (SELECT fp FROM knn_query)) AS dist
FROM knn_fps
ORDER BY dist, id
LIMIT 10;

-- the same search through the GiST index
SET enable_indexscan=on;
SET enable_seqscan=off;
CREATE TABLE knn_index AS
SELECT id, fp <%> (SELECT fp FROM knn_query) AS dist
FROM knn_fps
ORDER BY fp <%> (SELECT fp FROM knn_query)
LIMIT 10;

SELECT id, round(dist::numeric, 4) AS dist FROM knn_index ORDER BY dist, id;

-- compare distances, not ids: molecules at equal distance may come back in any order
SELECT (SELECT array_agg(round(dist::numeric, 6) ORDER BY dist) FROM knn_exact) =
       (SELECT array_agg(round(dist::numeric, 6) ORDER BY dist) FROM knn_index) AS same_top10;

DROP TABLE knn_exact, knn_index;


-- Dice

SET enable_indexscan=off;
SET enable_seqscan=on;
CREATE TABLE knn_exact AS
SELECT id, 1.0 - dice_sml(fp, (SELECT fp FROM knn_query)) AS dist
FROM knn_fps
ORDER BY dist, id
LIMIT 10;

SET enable_indexscan=on;
SET enable_seqscan=off;
CREATE TABLE knn_index AS
SELECT id, fp <#> (SELECT fp FROM knn_query) AS dist
FROM knn_fps
ORDER BY fp <#> (SELECT fp FROM knn_query)
LIMIT 10;

SELECT id, round(dist::numeric, 4) AS dist FROM knn_index ORDER BY dist, id;

SELECT (SELECT array_agg(round(dist::numeric, 6) ORDER BY dist) FROM knn_exact) =
       (SELECT array_agg(round(dist::numeric, 6) ORDER BY dist) FROM knn_index) AS same_top10;

RESET enable_indexscan;
RESET enable_bitmapscan;
RESET enable_seqscan;
RESET rdkit.morgan_fp_size;
DROP TABLE knn_mols, knn_fps, knn_query, knn_exact, knn_index;
