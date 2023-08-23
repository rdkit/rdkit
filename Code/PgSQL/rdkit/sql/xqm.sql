-- parsing and conversion to text
select mol_to_xqmol('Cc1n[nH]c(F)c1'::mol);
select mol_to_xqmol('COC |LN:1:1.3|'::mol);
select mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol);
select mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol, false);
select mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol, true, false);
select mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol, false, false);
select mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol, false, false, false);

-- substructure searching
select 'Cc1[nH]nc(F)c1'::mol @> mol_to_xqmol('Cc1n[nH]c(F)c1'::mol);
select substruct('Cc1[nH]nc(F)c1'::mol,mol_to_xqmol('Cc1n[nH]c(F)c1'::mol));
select substruct('Cc1[nH]nc(F)c1'::mol,mol_to_xqmol('Cc1n[nH]c(F)n1'::mol));

select substruct('COC'::mol,'COC |LN:1:1.3|'::mol);
select substruct('COC'::mol,mol_to_xqmol('COC |LN:1:1.3|'::mol));
select substruct('COOC'::mol,'COC |LN:1:1.3|'::mol);
select substruct('COOC'::mol,mol_to_xqmol('COC |LN:1:1.3|'::mol));
select substruct('COOOC'::mol,mol_to_xqmol('COC |LN:1:1.3|'::mol));
select substruct('COOOOC'::mol,mol_to_xqmol('COC |LN:1:1.3|'::mol));

select 'COc1[nH]ncc1'@>mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol);
select 'COc1[n][nH]cc1'@>mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol);
select 'COOc1[nH]ncc1'@>mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol);
select 'COOOc1[nH]ncc1'@>mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol);
select 'COOOOc1[nH]ncc1'@>mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol);


-- edge cases and error handling
select mol_to_xqmol('C1'::mol);
select mol_to_xqmol('c1cc1'::mol);
select mol_to_xqmol('COC1=NNC=C1'::mol);
select 'CCOC1=NNC=C1'::mol@>mol_to_xqmol('COC1=NNC=C1'::mol);
select mol_to_xqmol('CCCC'::mol);
select 'CCCCOC'::mol@>mol_to_xqmol('CCCC'::mol);

-- interaction with the generic groups
select 'COC1=NNC(CC)=C1'::mol @>> mol_adjust_query_properties('COC1=NNC(*)=C1 |$;;;;;;AEL_p;$,LN:1:1.3|'::mol,
  '{"makeDummiesQueries":true, "setGenericQueryFromProperties": true}');
select 'COC1=NNC(C=C)=C1'::mol @>> mol_adjust_query_properties('COC1=NNC(*)=C1 |$;;;;;;AEL_p;$,LN:1:1.3|'::mol,
  '{"makeDummiesQueries":true, "setGenericQueryFromProperties": true}');
select 'COC1=NNC(CC)=C1'::mol @>> mol_to_xqmol(mol_adjust_query_properties('COC1=NNC(*)=C1 |$;;;;;;AEL_p;$,LN:1:1.3|'::mol,
  '{"makeDummiesQueries":true, "setGenericQueryFromProperties": true}'));
select 'COC1=NNC(C=C)=C1'::mol @>> mol_to_xqmol(mol_adjust_query_properties('COC1=NNC(*)=C1 |$;;;;;;AEL_p;$,LN:1:1.3|'::mol,
  '{"makeDummiesQueries":true, "setGenericQueryFromProperties": true}'));
select 'COOC1=NNC(CC)=C1'::mol @>> mol_to_xqmol(mol_adjust_query_properties('COC1=NNC(*)=C1 |$;;;;;;AEL_p;$,LN:1:1.3|'::mol,
  '{"makeDummiesQueries":true, "setGenericQueryFromProperties": true}'));
select 'COOC1=NNC(C=C)=C1'::mol @>> mol_to_xqmol(mol_adjust_query_properties('COC1=NNC(*)=C1 |$;;;;;;AEL_p;$,LN:1:1.3|'::mol,
  '{"makeDummiesQueries":true, "setGenericQueryFromProperties": true}'));
