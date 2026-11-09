-- Binary COPY (send/recv) round-trip tests for mol, qmol, xqmol, and reaction

-- mol (mol_send / mol_recv)
CREATE TEMP TABLE binary_io_mol_src (id int, m mol);
INSERT INTO binary_io_mol_src VALUES
  (1, 'c1ccccc1'::mol),
  (2, 'C[C@@H]([C@H](C)F)[C@@H](C)Br'::mol),
  (3, NULL);
CREATE TEMP TABLE binary_io_mol_dst (id int, m mol);
\copy binary_io_mol_src TO 'results/binary_io_mol.bin' WITH (FORMAT binary)
\copy binary_io_mol_dst FROM 'results/binary_io_mol.bin' WITH (FORMAT binary)
SELECT d.id, d.m, (s.m @= d.m) AS same_mol
FROM binary_io_mol_src s INNER JOIN binary_io_mol_dst d USING (id)
ORDER BY d.id;
DROP TABLE binary_io_mol_src, binary_io_mol_dst;

-- qmol (qmol_send / qmol_recv)
CREATE TEMP TABLE binary_io_qmol_src (id int, q qmol);
INSERT INTO binary_io_qmol_src VALUES
  (1, 'c1cccc[n,c]1'::qmol),
  (2, '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'::qmol),
  (3, NULL);
CREATE TEMP TABLE binary_io_qmol_dst (id int, q qmol);
\copy binary_io_qmol_src TO 'results/binary_io_qmol.bin' WITH (FORMAT binary)
\copy binary_io_qmol_dst FROM 'results/binary_io_qmol.bin' WITH (FORMAT binary)
SELECT d.id, mol_to_smarts(d.q) AS smarts,
       (mol_to_smarts(s.q)::text = mol_to_smarts(d.q)::text) AS same_smarts,
       ('c1cccnc1'::mol @> d.q) AS matches_pyr,
       ('c1ccccc1'::mol @> d.q) AS matches_bz
FROM binary_io_qmol_src s INNER JOIN binary_io_qmol_dst d USING (id)
ORDER BY d.id;
DROP TABLE binary_io_qmol_src, binary_io_qmol_dst;

-- xqmol (xqmol_send / xqmol_recv)
CREATE TEMP TABLE binary_io_xqmol_src (id int, x xqmol);
INSERT INTO binary_io_xqmol_src VALUES
  (1, mol_to_xqmol('Cc1n[nH]c(F)c1'::mol)),
  (2, mol_to_xqmol('COC |LN:1:1.3|'::mol)),
  (3, NULL);
CREATE TEMP TABLE binary_io_xqmol_dst (id int, x xqmol);
\copy binary_io_xqmol_src TO 'results/binary_io_xqmol.bin' WITH (FORMAT binary)
\copy binary_io_xqmol_dst FROM 'results/binary_io_xqmol.bin' WITH (FORMAT binary)
SELECT d.id,
       (s.x::text = d.x::text) AS exact_match,
       ('Cc1[nH]nc(F)c1'::mol @> d.x) AS match_tautomer,
       ('COOC'::mol @> d.x) AS match_linknode
FROM binary_io_xqmol_src s INNER JOIN binary_io_xqmol_dst d USING (id)
ORDER BY d.id;
DROP TABLE binary_io_xqmol_src, binary_io_xqmol_dst;

-- reaction (reaction_send / reaction_recv)
CREATE TEMP TABLE binary_io_rxn_src (id int, r reaction);
INSERT INTO binary_io_rxn_src VALUES
  (1, reaction_from_smiles('c1ccccc1>>c1ccncc1')),
  (2, reaction_from_smiles('c1ccccc1>CC(=O)O>c1ccncc1')),
  (3, NULL);
CREATE TEMP TABLE binary_io_rxn_dst (id int, r reaction);
\copy binary_io_rxn_src TO 'results/binary_io_rxn.bin' WITH (FORMAT binary)
\copy binary_io_rxn_dst FROM 'results/binary_io_rxn.bin' WITH (FORMAT binary)
SELECT d.id, d.r, (s.r @= d.r) AS same_rxn
FROM binary_io_rxn_src s INNER JOIN binary_io_rxn_dst d USING (id)
ORDER BY d.id;
DROP TABLE binary_io_rxn_src, binary_io_rxn_dst;

-- malformed binary payloads raise ERRCODE_DATA_EXCEPTION cleanly
\copy (SELECT '\xdeadbeef'::bytea) TO 'results/binary_io_bad.bin' WITH (FORMAT binary)
CREATE TEMP TABLE binary_io_bad_mol (m mol);
\copy binary_io_bad_mol FROM 'results/binary_io_bad.bin' WITH (FORMAT binary)
DROP TABLE binary_io_bad_mol;

CREATE TEMP TABLE binary_io_bad_qmol (q qmol);
\copy binary_io_bad_qmol FROM 'results/binary_io_bad.bin' WITH (FORMAT binary)
DROP TABLE binary_io_bad_qmol;

CREATE TEMP TABLE binary_io_bad_xqmol (x xqmol);
\copy binary_io_bad_xqmol FROM 'results/binary_io_bad.bin' WITH (FORMAT binary)
DROP TABLE binary_io_bad_xqmol;

CREATE TEMP TABLE binary_io_bad_rxn (r reaction);
\copy binary_io_bad_rxn FROM 'results/binary_io_bad.bin' WITH (FORMAT binary)
DROP TABLE binary_io_bad_rxn;
