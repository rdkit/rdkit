select regexp_replace(mol_to_json('CC'::mol)::text,'"toolkitVersion":".*?"','"toolkitVersion":""');
select regexp_replace(mol_to_json('C[O,N]'::qmol)::text,'"toolkitVersion":".*?"','"toolkitVersion":""');
