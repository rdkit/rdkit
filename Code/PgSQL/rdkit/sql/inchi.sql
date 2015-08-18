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

