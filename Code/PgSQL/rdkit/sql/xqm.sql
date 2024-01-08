-- parsing and conversion to text
select mol_to_xqmol('Cc1n[nH]c(F)c1'::mol);
select mol_to_xqmol('COC |LN:1:1.3|'::mol);
select mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol);
select mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol, false);
select mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol, true, false);
select mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol, false, false);
select mol_to_xqmol('COC1=NNC=C1 |LN:1:1.3|'::mol, false, false, true);

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

-- problems we've encountered
select 'CCOC(=O)c1cc2n(n1)C(C(=O)O)Nc1cc(Cl)ccc1-2'::mol @> mol_to_xqmol('COCc1n[nH]c(C)c1 |LN:1:1.3|'::mol);

-- link nodes and variable attachment points
select 'Cc1ccnc2nc(CN)[nH]c12'::mol@>>mol_to_xqmol(mol_from_ctab('qry 
  Mrv2305 09052314502D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 13 13 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -4.75 1.9567 0 0
M  V30 2 C -6.0837 1.1867 0 0
M  V30 3 C -6.0837 -0.3534 0 0
M  V30 4 C -4.75 -1.1234 0 0
M  V30 5 C -3.4163 -0.3534 0 0
M  V30 6 C -3.4163 1.1867 0 0
M  V30 7 N -1.9692 1.7134 0 0
M  V30 8 N -1.8822 -0.7768 0 0
M  V30 9 C -1.0211 0.4999 0 0
M  V30 10 C 0.5179 0.5536 0 0
M  V30 11 N 1.2409 1.9133 0 0
M  V30 12 * -5.6391 -0.0967 0 0
M  V30 13 C -5.6391 -2.4067 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 8 9
M  V30 8 1 7 6
M  V30 9 1 5 8
M  V30 10 2 7 9
M  V30 11 1 9 10
M  V30 12 1 10 11
M  V30 13 1 12 13 ENDPTS=(3 4 3 2) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END'));
select 'c1ccnc2nc(CN)[nH]c12'::mol@>>mol_to_xqmol(mol_from_ctab('qry 
  Mrv2305 09052314502D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 13 13 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -4.75 1.9567 0 0
M  V30 2 C -6.0837 1.1867 0 0
M  V30 3 C -6.0837 -0.3534 0 0
M  V30 4 C -4.75 -1.1234 0 0
M  V30 5 C -3.4163 -0.3534 0 0
M  V30 6 C -3.4163 1.1867 0 0
M  V30 7 N -1.9692 1.7134 0 0
M  V30 8 N -1.8822 -0.7768 0 0
M  V30 9 C -1.0211 0.4999 0 0
M  V30 10 C 0.5179 0.5536 0 0
M  V30 11 N 1.2409 1.9133 0 0
M  V30 12 * -5.6391 -0.0967 0 0
M  V30 13 C -5.6391 -2.4067 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 8 9
M  V30 8 1 7 6
M  V30 9 1 5 8
M  V30 10 2 7 9
M  V30 11 1 9 10
M  V30 12 1 10 11
M  V30 13 1 12 13 ENDPTS=(3 4 3 2) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END'));