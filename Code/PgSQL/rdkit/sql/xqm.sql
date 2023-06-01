-- parsing and conversion to text
select mol_to_tautomerquery('Cc1n[nH]c(F)c1'::mol);
select mol_enumeratequery('COC |LN:1:1.3|'::mol);
select xqmol_to_tautomerquery(mol_enumeratequery('COC1=NNC=C1 |LN:1:1.3|'::mol));

-- substructure searching
select 'Cc1[nH]nc(F)c1'::mol @> mol_to_tautomerquery('Cc1n[nH]c(F)c1'::mol);
select substruct('Cc1[nH]nc(F)c1'::mol,mol_to_tautomerquery('Cc1n[nH]c(F)c1'::mol));
select substruct('Cc1[nH]nc(F)c1'::mol,mol_to_tautomerquery('Cc1n[nH]c(F)n1'::mol));

select substruct('COC'::mol,'COC |LN:1:1.3|'::mol);
select substruct('COC'::mol,mol_enumeratequery('COC |LN:1:1.3|'::mol));
select substruct('COOC'::mol,'COC |LN:1:1.3|'::mol);
select substruct('COOC'::mol,mol_enumeratequery('COC |LN:1:1.3|'::mol));
select substruct('COOOC'::mol,mol_enumeratequery('COC |LN:1:1.3|'::mol));
select substruct('COOOOC'::mol,mol_enumeratequery('COC |LN:1:1.3|'::mol));

select 'COc1[nH]ncc1'@>xqmol_to_tautomerquery(mol_enumeratequery('COC1=NNC=C1 |LN:1:1.3|'::mol));
select 'COc1[n][nH]cc1'@>xqmol_to_tautomerquery(mol_enumeratequery('COC1=NNC=C1 |LN:1:1.3|'::mol));
select 'COOc1[nH]ncc1'@>xqmol_to_tautomerquery(mol_enumeratequery('COC1=NNC=C1 |LN:1:1.3|'::mol));
select 'COOOc1[nH]ncc1'@>xqmol_to_tautomerquery(mol_enumeratequery('COC1=NNC=C1 |LN:1:1.3|'::mol));
select 'COOOOc1[nH]ncc1'@>xqmol_to_tautomerquery(mol_enumeratequery('COC1=NNC=C1 |LN:1:1.3|'::mol));


-- edge cases and error handling
select mol_to_tautomerquery('C1'::mol);
select mol_to_tautomerquery('c1cc1'::mol);
select mol_enumeratequery('COC1=NNC=C1'::mol);
select 'CCOC1=NNC=C1'::mol@>mol_enumeratequery('COC1=NNC=C1'::mol);
select mol_to_tautomerquery('CCCC'::mol);
select 'CCCCOC'::mol@>mol_to_tautomerquery('CCCC'::mol);

