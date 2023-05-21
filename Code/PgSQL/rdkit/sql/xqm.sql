select 'Cc1[nH]nc(F)c1'::mol @> mol_to_tautomerquery('Cc1n[nH]c(F)c1'::mol);
select substruct('Cc1[nH]nc(F)c1'::mol,mol_to_tautomerquery('Cc1n[nH]c(F)c1'::mol));
select substruct('Cc1[nH]nc(F)c1'::mol,mol_to_tautomerquery('Cc1n[nH]c(F)n1'::mol));

select substruct('COC'::mol,'COC |LN:1:1.3|'::mol);
select substruct('COC'::mol,mol_enumeratequery('COC |LN:1:1.3|'::mol));
select substruct('COOC'::mol,'COC |LN:1:1.3|'::mol);
select substruct('COOC'::mol,mol_enumeratequery('COC |LN:1:1.3|'::mol));

-- error handling
select mol_to_tautomerquery('C1'::mol);
select mol_to_tautomerquery('c1cc1'::mol);
