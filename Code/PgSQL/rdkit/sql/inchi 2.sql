select mol_inchi('c1ccccc1'::mol);
select mol_inchikey('c1ccccc1'::mol);
select mol_inchi('Cc1cc(C)[n+]c(C)c1'::mol);
select mol_inchikey('Cc1cc(C)[n+]c(C)c1'::mol);
select mol_inchi(mol_from_ctab((Chr(10) || Chr(10) || Chr(10) ||
'  0  0  0  0  0  0  0  0  0  0999 V2000' ||
Chr(10) ||
'M  END')::cstring));
select mol_inchi('');
select mol_inchikey('');
select mol_inchi('CC1=NC=CN1'::mol);
select mol_inchikey('CC1=NC=CN1'::mol);
select mol_inchi('CC1=NC=CN1'::mol,'/FixedH');
select mol_inchikey('CC1=NC=CN1'::mol,'/FixedH');

-- Non-InChI-able molecules should return NULL.
select coalesce(mol_inchi('CC*'), '<NULL>');
select coalesce(mol_inchikey('CC*'), '<NULL>');
