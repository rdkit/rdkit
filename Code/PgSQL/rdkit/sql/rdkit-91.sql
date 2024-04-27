--
-- first, define the datatype.  Turn off echoing so that expected file
-- does not depend on contents of rdkit.sql.
--
SET client_min_messages = warning;
\set ECHO none
CREATE EXTENSION rdkit;
\set ECHO all
RESET client_min_messages;
SET extra_float_digits=0;

SELECT is_valid_smiles('c1ccccc1');
SELECT mol_from_smiles('c1ccccc1');
SELECT is_valid_smiles('c1ccccc');
SELECT mol_from_smiles('c1ccccc');
SELECT mol_from_smiles('cccccc');
SELECT is_valid_smiles('c1cccn1');
SELECT is_valid_smarts('c1ccc[n,c]1');
SELECT qmol_from_smarts('c1ccc[n,c]1');
SELECT is_valid_smarts('c1ccc');
SELECT qmol_from_smarts('c1ccc');
SELECT mol_to_smiles(mol_from_smiles('c1ccccc1'));
SELECT mol_to_smarts(mol_from_smiles('c1ccccc1'));
SELECT mol_to_smarts('c1cccc[n,c]1'::qmol);
SELECT mol_to_smiles('c1cccc[n,c]1'::qmol);
SELECT is_valid_smiles('');
SELECT mol_from_smiles('');
SELECT mol_to_smiles(mol_from_smiles(''));
SELECT mol_to_smiles(mol_from_smiles('C[C@H](F)[C@H](C)[C@@H](C)Br'));
SELECT mol_to_smiles(mol_from_smiles('C[C@H](F)[C@H](C)[C@@H](C)Br'), false);

CREATE TABLE pgmol (id int, m mol);
\copy pgmol from 'data/data'

CREATE UNIQUE INDEX mol_ididx ON pgmol (id);

SELECT count(*) FROM pgmol;

SELECT count(*) FROM pgmol WHERE m @> 'c1ccccc1';
SELECT count(*) FROM pgmol WHERE m @> 'c1cccnc1';
SELECT count(*) FROM pgmol WHERE 'c1ccccc1' <@ m;
SELECT count(*) FROM pgmol WHERE 'c1cccnc1' <@ m;

SELECT count(*) FROM pgmol WHERE m @> qmol_from_smarts('c1ccccc1');
SELECT count(*) FROM pgmol WHERE m @> qmol_from_smarts('c1cccnc1');
SELECT count(*) FROM pgmol WHERE m @> qmol_from_smarts('c1ccc[n,c]c1');
SELECT count(*) FROM pgmol WHERE qmol_from_smarts('c1ccccc1') <@ m;
SELECT count(*) FROM pgmol WHERE qmol_from_smarts('c1ccc[n,c]c1') <@ m;


SELECT id, rdkit_fp(m) AS f, maccs_fp(m) as maccsf INTO pgbfp FROM pgmol;
CREATE UNIQUE INDEX bfp_ididx ON pgbfp (id);

SELECT id, morgan_fp(m,1) AS f INTO pgsfp FROM pgmol;
CREATE UNIQUE INDEX sfp_ididx ON pgsfp (id);

SELECT id, torsion_fp(m) AS f INTO pgtorsfp FROM pgmol;
SELECT id, atompair_fp(m) AS f INTO pgpairfp FROM pgmol;

set rdkit.tanimoto_threshold=0.5;
set rdkit.dice_threshold=0.5;

SELECT
	id,
	tanimoto_sml(rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol), f)
FROM
	 (SELECT * FROM pgbfp ORDER BY id) AS t
WHERE rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol) % f
LIMIT 10;

SELECT
	id,
	dice_sml(rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol), f)
FROM
	 (SELECT * FROM pgbfp ORDER BY id) AS t
WHERE rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol) % f
LIMIT 10;

SELECT
	id,
	tanimoto_sml(rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol), f)
FROM
	 (SELECT * FROM pgbfp ORDER BY id) AS t
WHERE rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol) # f
LIMIT 10;

SELECT
	id,
	dice_sml(rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol), f),
	size(f)
FROM
	 (SELECT * FROM pgbfp ORDER BY id) AS t
WHERE rdkit_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol) # f
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

select dice_sml(morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol, 1), morgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)N)'::mol, 1)) sml;

select dice_sml(featmorgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol, 1), featmorgan_fp('C1C(OC2=CC(=CC(=C2C1=O)O)N)'::mol, 1)) sml;

select dice_sml(morganbv_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol, 1), morganbv_fp('C1C(OC2=CC(=CC(=C2C1=O)O)N)'::mol, 1)) sml;

select dice_sml(featmorganbv_fp('C1C(OC2=CC(=CC(=C2C1=O)O)O)'::mol, 1), featmorganbv_fp('C1C(OC2=CC(=CC(=C2C1=O)O)N)'::mol, 1)) sml;

select 'Cc1ccccc1'::mol@='c1ccccc1C'::mol;
select 'Cc1ccccc1'::mol@='c1ccccc1CC'::mol;
select 'Cc1ccccc1'::mol@='c1cccnc1C'::mol;

select subtract(torsion_fp('CCC1CCNCC1'),torsion_fp('OCC1CCNCC1'))=subtract(torsion_fp('CCC1CCOCC1'),torsion_fp('OCC1CCOCC1'));
select subtract(torsion_fp('CCC1CCNCC1'),torsion_fp('OCC1CCNCC1'))=subtract(torsion_fp('CCC1CCOCC1'),torsion_fp('NCC1CCOCC1'));
select add(torsion_fp('CCC1CCNCC1'),torsion_fp('OCC1CCNCC1'))=add(torsion_fp('CCC1CCOCC1'),torsion_fp('OCC1CCOCC1'));
select add(torsion_fp('CCC1CCNCC1'),torsion_fp('OCC1CCNCC1'))=add(torsion_fp('CCC1CCOCC1'),torsion_fp('NCC1CCOCC1'));

select add(torsion_fp('CCC1CCNCC1'),torsion_fp('OCC1CCNCC1'))=subtract(torsion_fp('CCC1CCNCC1'),torsion_fp('OCC1CCNCC1'));
select add(torsion_fp('CCC1CCNCC1'),torsion_fp('OCC1CCNCC1'))=subtract(torsion_fp('CCC1CCOCC1'),torsion_fp('OCC1CCOCC1'));

select is_valid_ctab('chiral1.mol
  ChemDraw04200416412D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -0.0141    0.0553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8109    0.0553    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4266    0.7697    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.0141   -0.7697    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.8109   -0.1583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  1
  1  5  1  0
M  END');
select is_valid_ctab('invalid');
select mol_from_ctab('chiral1.mol
  ChemDraw04200416412D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -0.0141    0.0553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8109    0.0553    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4266    0.7697    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.0141   -0.7697    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.8109   -0.1583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  1
  1  5  1  0
M  END');

-- mol_to_ctab() - suppress auto-generation of depiction.
select mol_to_ctab(mol('CCC'), false);
-- mol_to_ctab() - with auto-generated depiction.
select mol_to_ctab(mol('CCC'));
-- mol_to_ctab() - should preserve existing/input depiction. Note the
-- extra 'true' parameter to 'mol_from_ctab()' that forces the cartridge
-- to preserve the input conformer. Otherwise the conformer will be lost.
select mol_to_ctab(mol_from_ctab('chiral1.mol
  ChemDraw04200416412D

  5  4  0  0  0  0  0  0  0  0999 V2000
   -0.0141    0.0553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8109    0.0553    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4266    0.7697    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.0141   -0.7697    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.8109   -0.1583    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  1
  1  5  1  0
M  END', true));
-- mol_to_ctab() - force v3000.
select mol_to_ctab(mol('CCC'), false, true);
select mol_to_v3kctab(mol('CCC'), false);

select all_values_lt(torsion_fp('c1ccccc1C'::mol),2);
select all_values_lt(torsion_fp('c1ccccc1C'::mol),3);
select all_values_gt(torsion_fp('c1ccccc1C'::mol),1);
select all_values_gt(torsion_fp('c1ccccc1C'::mol),2);

select is_valid_mol_pkl('foo'::bytea);
select is_valid_mol_pkl(mol_to_pkl('c1ccccc1'::mol));
select mol_from_pkl(mol_to_pkl('c1ccccc1'::mol));

select tanimoto_sml(morganbv_fp('c1ccccn1'::mol),morganbv_fp('c1ccccc1'::mol));
select tanimoto_sml(bfp_from_binary_text(bfp_to_binary_text(morganbv_fp('c1ccccn1'::mol))),
                    bfp_from_binary_text(bfp_to_binary_text(morganbv_fp('c1ccccc1'::mol))));

-- GitHub issue 9
select 'C1CC2CC3C45C2C2C6C7C8C9C%10C(C1)C1C%11%10C%109C98C87C76C42C24C65C3C3C56C64C4%12C72C28C79C8%10C9%11C1C1C%109C98C87C42C24C7%12C%116C65C3C3C56C6%11C%117C74C4%12C82C29C8%10C1C1C98C42C24C89C1C1C98C84C4%10C%122C27C7%11C%116C65C3C3C56C6%11C%117C42C24C7%11C%116C65C3C3C56C6%11C%117C74C4%12C%102C28C89C1C1C98C42C24C89C1C1C98C84C4%10C%122C27C7%11C%116C65C3C3C56C6%11C%117C42C24C7%11C%116C65C3C3C56C6%11C%117C74C4%12C%102C28C89C1C1C98C42C24C89C1CC8C4C1C%122C27C4%11C76C65C3CC6C7C4C12'::mol;

-- chiral matching
set rdkit.do_chiral_sss=false;
select 'C[C@H](F)Cl'::mol@>'CC(F)Cl'::mol as match;
select 'C[C@H](F)Cl'::mol@>'C[C@H](F)Cl'::mol as match;
select 'C[C@H](F)Cl'::mol@>'C[C@@H](F)Cl'::mol as match;
set rdkit.do_chiral_sss=true;
select 'C[C@H](F)Cl'::mol@>'CC(F)Cl'::mol as match;
select 'C[C@H](F)Cl'::mol@>'C[C@H](F)Cl'::mol as match;
select 'C[C@H](F)Cl'::mol@>'C[C@@H](F)Cl'::mol as match;
set rdkit.do_chiral_sss=false;

-- github #2790
set rdkit.do_chiral_sss=false;
select 'C[C@H](F)Cl'::mol@='C[C@H](F)Cl'::mol;
select 'C[C@H](F)Cl'::mol@='C[C@@H](F)Cl'::mol;
select 'C[C@H](F)Cl'::mol@='CC(F)Cl'::mol;
select 'CC(F)Cl'::mol@='C[C@@H](F)Cl'::mol;
select 'CC(F)Cl'::mol@='C[C@H](F)Cl'::mol;
select 'CC(F)Cl'::mol@='CC(F)Cl'::mol;
set rdkit.do_chiral_sss=true;
select 'C[C@H](F)Cl'::mol@='C[C@H](F)Cl'::mol;
select 'C[C@H](F)Cl'::mol@='C[C@@H](F)Cl'::mol;
select 'C[C@H](F)Cl'::mol@='CC(F)Cl'::mol;
select 'CC(F)Cl'::mol@='C[C@@H](F)Cl'::mol;
select 'CC(F)Cl'::mol@='C[C@H](F)Cl'::mol;
select 'CC(F)Cl'::mol@='CC(F)Cl'::mol;
set rdkit.do_chiral_sss=false;

-- Enhanced stereo
set rdkit.do_chiral_sss=false;
set rdkit.do_enhanced_stereo_sss=false; /* has no effect when do_chiral_sss is false */
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol;
select 'C[C@H](O)[C@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol;
select 'C[C@H](O)[C@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@H](O)[C@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@H](O)[C@H](C)F'::mol;
set rdkit.do_enhanced_stereo_sss=true; /* has no effect when do_chiral_sss is false */
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol;
select 'C[C@H](O)[C@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol;
select 'C[C@H](O)[C@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@H](O)[C@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@H](O)[C@H](C)F'::mol;
set rdkit.do_chiral_sss=true;
set rdkit.do_enhanced_stereo_sss=false;
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol;
select 'C[C@H](O)[C@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol;
select 'C[C@H](O)[C@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@H](O)[C@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@H](O)[C@H](C)F'::mol;
set rdkit.do_enhanced_stereo_sss=true; /* now we expect to see an effect */
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol;
select 'C[C@H](O)[C@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol;
select 'C[C@H](O)[C@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@@H](O)[C@@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol@>'C[C@H](O)[C@H](C)F'::mol;
select 'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol@>'C[C@H](O)[C@H](C)F'::mol;

set rdkit.do_chiral_sss=false;
set rdkit.do_enhanced_stereo_sss=false;

-- forcing chiral queries
select substruct_chiral('C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol,'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol);
select substruct_chiral('C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol,'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol);
select rsubstruct_chiral('C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol,'C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol);
select rsubstruct_chiral('C[C@@H](O)[C@@H](C)F |o1:1,3,r|'::mol,'C[C@@H](O)[C@@H](C)F |&1:1,3,r|'::mol);


-- substructure counts
select substruct_count('c1ccncc1'::mol,'c1ccncc1'::mol);
select substruct_count('c1ccncc1'::mol,'c1ccncc1'::mol,false);
select substruct_count('c1ccccc1C[C@@H](O)[C@@H](C)F |&1:7,9,r|'::mol,'c1ccccc1C[C@@H](O)[C@@H](C)F |o1:7,9,r|'::mol);
select substruct_count('c1ccccc1C[C@@H](O)[C@@H](C)F |o1:7,9,r|'::mol,'c1ccccc1C[C@@H](O)[C@@H](C)F |&1:7,9,r|'::mol);
select substruct_count_chiral('c1ccccc1C[C@@H](O)[C@@H](C)F |&1:7,9,r|'::mol,'c1ccccc1C[C@@H](O)[C@@H](C)F |o1:7,9,r|'::mol);
select substruct_count_chiral('c1ccccc1C[C@@H](O)[C@@H](C)F |o1:7,9,r|'::mol,'c1ccccc1C[C@@H](O)[C@@H](C)F |&1:7,9,r|'::mol);
select substruct_count_chiral('c1ccccc1C[C@@H](O)[C@@H](C)F |&1:7,9,r|'::mol,'c1ccccc1C[C@@H](O)[C@@H](C)F |o1:7,9,r|'::mol,false);
select substruct_count_chiral('c1ccccc1C[C@@H](O)[C@@H](C)F |o1:7,9,r|'::mol,'c1ccccc1C[C@@H](O)[C@@H](C)F |&1:7,9,r|'::mol,false);

-- special queries
select 'c1ccc[nH]1'::mol@>mol_from_smiles('c1cccn1[H]') as match;
select 'c1cccn1C'::mol@>mol_from_smiles('c1cccn1[H]') as match;
select 'c1ccc[nH]1'::mol@>qmol_from_smiles('c1cccn1[H]') as match;
select 'c1cccn1C'::mol@>qmol_from_smiles('c1cccn1[H]') as match;
select 'c1ccc[nH]1'::mol@>mol_from_ctab('query
  Mrv0541 04021509592D

  6  6  0  0  0  0            999 V2000
   -0.2652    0.7248    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2652   -0.1002    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  4  5  2  0  0  0  0
  1  5  1  0  0  0  0
  3  4  1  0  0  0  0
  1  6  1  0  0  0  0
M  END') as match;
select 'c1cccn1C'::mol@>mol_from_ctab('query
  Mrv0541 04021509592D

  6  6  0  0  0  0            999 V2000
   -0.2652    0.7248    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2652   -0.1002    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  4  5  2  0  0  0  0
  1  5  1  0  0  0  0
  3  4  1  0  0  0  0
  1  6  1  0  0  0  0
M  END') as match;
select 'c1ccc[nH]1'::mol@>qmol_from_ctab('query
  Mrv0541 04021509592D

  6  6  0  0  0  0            999 V2000
   -0.2652    0.7248    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2652   -0.1002    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  4  5  2  0  0  0  0
  1  5  1  0  0  0  0
  3  4  1  0  0  0  0
  1  6  1  0  0  0  0
M  END') as match;
select 'c1cccn1C'::mol@>qmol_from_ctab('query
  Mrv0541 04021509592D

  6  6  0  0  0  0            999 V2000
   -0.2652    0.7248    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2652   -0.1002    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  4  5  2  0  0  0  0
  1  5  1  0  0  0  0
  3  4  1  0  0  0  0
  1  6  1  0  0  0  0
M  END') as match;
-- github #4787:
select mol_to_smarts(qmol_from_ctab('query
  Mrv0541 04021509592D

  6  6  0  0  0  0            999 V2000
   -0.2652    0.7248    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2652   -0.1002    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  4  5  2  0  0  0  0
  1  5  1  0  0  0  0
  3  4  1  0  0  0  0
  1  6  1  0  0  0  0
M  END'));
select mol_to_smarts(qmol_from_ctab('query
  Mrv0541 04021509592D

  6  6  0  0  0  0            999 V2000
   -0.2652    0.7248    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9796    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.9623    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4493    1.1373    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2652   -0.1002    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  4  5  4  0  0  0  0
  1  5  4  0  0  0  0
  3  4  4  0  0  0  0
  1  6  1  0  0  0  0
M  END'));
select mol_to_smarts(qmol_from_ctab('Boronate acid/ester(aryl)
  SciTegic12012112112D

  5  4  0  0  0  0            999 V2000
    1.7243   -2.7324    0.0000 A   0  0
    2.7559   -2.1456    0.0000 C   0  0
    3.7808   -2.7324    0.0000 B   0  0
    4.8057   -2.1456    0.0000 O   0  0
    3.7808   -3.9190    0.0000 O   0  0
  1  2  4  0  0  1  0
  2  3  1  0
  3  4  1  0
  3  5  1  0
M  END'));
select mol_to_smarts(qmol_from_smiles('c:c'));
select mol_to_smarts(qmol_from_smiles('C1=CC=CC=C1'));

-- mol_adjust_query_properties
select 'C1CC1OC'::mol @> 'C1CC1O*'::mol;
select 'C1CC1OC'::mol @> mol_adjust_query_properties('C1CC1O*'::mol);
select 'C1CC1OC'::mol @> mol_adjust_query_properties('C1CC1O*'::mol,'{"makeDummiesQueries":false}');
select 'CC1CC1OC'::mol @> 'C1CC1O*'::mol;
select 'CC1CC1OC'::mol @> mol_adjust_query_properties('C1CC1O*'::mol);
select 'CC1CC1OC'::mol @> mol_adjust_query_properties('C1CC1O*'::mol,'{"adjustDegree":false}');
select 'C1CC1C(C)OC'::mol @> mol_adjust_query_properties('C1CC1CO*'::mol);
select 'C1CC1C(C)OC'::mol @> mol_adjust_query_properties('C1CC1CO*'::mol,'{"adjustDegreeFlags":"IGNOREDUMMIES"}');
select 'C1CC1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol);
select 'C1CC1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"IGNORENONE"}');
select 'C1CC1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"IGNOREALL"}');
select 'C1CC1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"IGNORECHAINS"}');
select 'C1CC1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"IGNORERINGS"}');
select 'C1CC1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"IGNORERINGS|IGNORECHAINS"}');
select 'C1C(C)C1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol);
select 'C1C(C)C1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"IGNORENONE"}');
select 'C1C(C)C1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"IGNOREALL"}');
select 'C1C(C)C1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"IGNORECHAINS"}');
select 'C1C(C)C1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"IGNORERINGS"}');
select 'C1C(C)C1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"IGNORERINGS|IGNORECHAINS"}');
select 'C1C(C)C1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegreeFlags":"bogus"}');
select 'C1C([2H])C1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol);
select 'C1C([2H])C1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegree":false}');
select 'C1C([2H])C1CCCC'::mol @> mol_adjust_query_properties('C1CC1CC'::mol,'{"adjustDegree":false,"adjustHeavyDegree":true}');
select mol_to_smarts(mol_adjust_query_properties('*c1ncc(*)cc1'::mol));
select mol_to_smarts(mol_adjust_query_properties('*c1ncc(*)cc1'::qmol));


-- CXSmiles
SELECT mol_to_smiles(mol_from_smiles('C[C@H](F)[C@H](C)[C@@H](C)Br |a:1,o1:3,5|'));
SELECT mol_to_cxsmiles(mol_from_smiles('C[C@H](F)[C@H](C)[C@@H](C)Br |a:1,o1:3,5|'));
SELECT mol_to_cxsmiles(mol_from_smiles('C[C@H](F)[C@H](C)[C@@H](C)Br |a:1,o1:3,5|'), false);
SELECT mol_to_cxsmarts(mol_from_smiles('C[C@H](F)[C@H](C)[C@@H](C)Br |a:1,o1:3,5|'));
SELECT mol_to_cxsmarts(qmol_from_smarts('C[C@H]([F,Cl,Br])[C@H](C)[C@@H](C)Br |a:1,o1:3,5|'));

-- CXSmiles from mol_out
SELECT mol_out(mol_from_smiles('C[C@H](F)[C@H](C)[C@@H](C)Br |a:1,o1:3,5|'));

-- github #3688: bad input to qmol_from_ctab() crashes db
select qmol_from_ctab('a'::cstring,false);
-- github #3689: bad input to qmol_from_smiles() crashes db
select qmol_from_smiles('a'::cstring);
select qmol_from_smiles('C1C'::cstring);

-- casting from mol to qmol
select mol_from_smiles('C=C')::qmol;

-- github #5095: cannot restore molecule
select mol_in('c1cccc'::cstring);
select mol_in('c1cccc1'::cstring);
select mol_in('c1co(C)cc1'::cstring);
select mol_in('c1cccc'::cstring);
select mol_in('CN(=O)=O'::cstring);
select 'CN(=O)=O'::mol;
select 'c1cccc1'::mol;
select 'c1co(C)cc1'::mol;
select mol_in('c1cccc1'::cstring) @> '[r5]'::qmol;
select 'c1cccc1'::mol @> '[r5]'::qmol;
select mol_in('Cc1ccc2c(c1)-n1-c(=O)c=cc(=O)-n-2-c2cc(C)ccc2-1');
select 'c1cccc1'::text::mol;
select 'c1cccc1'::varchar::mol;
select mol_from_smiles('CCN(=O)=O') @> 'CN(=O)=O';

-- github #6002: molcmp failure
select mol_cmp(mol_from_ctab('
  Mrv2211 02092314292D

  5  4  0  0  0  0            999 V2000
    0.0000    3.6020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    4.0145    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    1.4290    4.4270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1270    3.3001    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3020    4.7291    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  6  0  0  0  0
  2  3  2  0  0  0  0
  2  4  2  0  0  0  0
  2  5  6  0  0  0  0
M  END'),mol_from_ctab('
  Mrv2211 02092314292D

  5  4  0  0  0  0            999 V2000
    0.0000    3.6020    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    4.0145    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    1.4290    4.4270    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.1270    3.3001    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.3020    4.7291    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  6  0  0  0  0
  2  3  2  0  0  0  0
  2  4  2  0  0  0  0
  2  5  6  0  0  0  0
M  END'));

-- mol properties being properly serialized
select 'COC1=NNC=C1 |LN:1:1.3|'::mol;

-- bond properties all preserved
select mol_to_v3kctab(mol_from_ctab('
  Mrv2211 09062306242D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -7.7917 4.0833 0 0
M  V30 2 C -6.458 4.8533 0 0
M  V30 3 C -5.1243 4.0833 0 0
M  V30 4 * -6.458 4.34 0 0
M  V30 5 C -5.303 6.3405 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 4 5 ENDPTS=(3 1 2 3) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END'));