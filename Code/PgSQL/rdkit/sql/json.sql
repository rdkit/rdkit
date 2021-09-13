select regexp_replace(mol_to_json('CC'::mol)::text,'"toolkitVersion":".*?"','"toolkitVersion":""');
select regexp_replace(mol_to_json('C[O,N]'::qmol)::text,'"toolkitVersion":".*?"','"toolkitVersion":""');
select mol_to_smiles(mol_from_json(mol_to_json('CC'::mol)));
select mol_to_smarts(mol_from_json(mol_to_json('C[O,N]'::qmol)));
select mol_to_smarts(qmol_from_json(mol_to_json('C[O,N]'::qmol)));
