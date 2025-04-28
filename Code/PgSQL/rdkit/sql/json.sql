select regexp_replace(mol_to_json('CC'::mol)::text,'ersion":.*?([,}])','ersion":""\1','g');
select regexp_replace(mol_to_json('C[O,N]'::qmol)::text,'ersion":.*?([,}])','ersion":""\1','g');
select mol_to_smiles(mol_from_json(mol_to_json('CC'::mol)));
select mol_to_smarts(mol_from_json(mol_to_json('C[O,N]'::qmol)));
select mol_to_smarts(qmol_from_json(mol_to_json('C[O,N]'::qmol)));
