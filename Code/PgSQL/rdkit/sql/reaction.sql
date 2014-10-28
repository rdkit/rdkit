SET rdkit.ignore_reaction_agents=false;
SET rdkit.agent_FP_bit_ratio=0.2;
SET rdkit.difference_FP_weight_agents=1;
SET rdkit.difference_FP_weight_nonagents=10;
SET rdkit.move_unmmapped_reactants_to_agents=true;
SET rdkit.threshold_unmapped_reactant_atoms=0.2;
SET rdkit.init_reaction=true;

SELECT reaction_from_smiles('c1ccccc1>>c1cccnc1');
SELECT reaction_from_smiles('c1ccccc1>CC(=O)O>c1cccnc1');
SELECT reaction_from_smarts('[c1:1][c:2][c:3][c:4]c[c1:5]>CC(=O)O>[c1:1][c:2][c:3][c:4]n[c1:5]');
SELECT reaction_from_smarts('C(F)(F)F.[c1:1][c:2][c:3][c:4]c[c1:5]>CC(=O)O>[c1:1][c:2][c:3][c:4]n[c1:5]');
SELECT reaction_from_smarts('c1ccc[n,c]c1>>c1nccnc1');
SELECT reaction_to_smiles(reaction_from_smiles('c1ccccc1>>c1cccnc1'));
SELECT reaction_to_smiles(reaction_from_smarts('c1ccc[n,c]c1>>c1nccnc1'));
SELECT reaction_to_smarts(reaction_from_smarts('c1ccc[n,c]c1>>c1nccnc1'));
SELECT reaction_to_smarts('c1cccnc1>>c1nccnc1'::reaction);
SELECT reaction_to_ctab(reaction_from_smiles('c1ccccc1>>c1cccnc1'));
SELECT reaction_numreactants(reaction_from_smiles('[Cl].c1ccccc1>>c1cccnc1.[OH2]'));
SELECT reaction_numproducts(reaction_from_smiles('[Cl].c1ccccc1>>c1cccnc1.[OH2]'));
SELECT reaction_numagents(reaction_from_smiles('[Cl].c1ccccc1>CC(=O)O.[Na+]>c1cccnc1.[OH2]'));
SELECT reaction_numagents(reaction_from_smarts('C(F)(F)F.[c1:1][c:2][c:3][c:4]c[c1:5]>CC(=O)O>[c1:1][c:2][c:3][c:4]n[c1:5]'));
SET rdkit.move_unmmapped_reactants_to_agents=false;
SELECT reaction_numagents(reaction_from_smarts('C(F)(F)F.[c1:1][c:2][c:3][c:4]c[c1:5]>CC(=O)O>[c1:1][c:2][c:3][c:4]n[c1:5]'));
SET rdkit.move_unmmapped_reactants_to_agents=true;
SET rdkit.threshold_unmapped_reactant_atoms=0.9;
SELECT reaction_numagents(reaction_from_smarts('C(F)(F)F.[c1:1][c:2][c:3][c:4]c[c1:5]>CC(=O)O>[c1:1][c:2][c:3][c:4]n[c1:5]'));
SET rdkit.threshold_unmapped_reactant_atoms=0.2;
SELECT 'c1ccccc1>>c1cccnc1'::reaction @= 'c1ccccc1>>c1cccnc1'::reaction;
SELECT 'c1ccccc1>>c1cccnc1'::reaction @= 'c1ccccc1>>c1cncnc1'::reaction;
SELECT reaction_from_ctab('$RXN                                                                 
                                                                     
      RDKit                                                          
                                                                     
  1  1                                                               
$MOL                                                                 
                                                                     
     RDKit                                                           
                                                                     
  6  6  0  0  0  0  0  0  0  0999 V2000                              
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0                                                         
  2  3  4  0                                                         
  3  4  4  0                                                         
  4  5  4  0                                                         
  5  6  4  0                                                         
  6  1  4  0                                                         
M  END                                                               
$MOL                                                                 
                                                                     
     RDKit                                                           
                                                                     
  6  6  0  0  0  0  0  0  0  0999 V2000                              
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0                                                         
  2  3  4  0                                                         
  3  4  4  0                                                         
  4  5  4  0                                                         
  5  6  4  0                                                         
  6  1  4  0                                                         
M  END');

CREATE TABLE tmp (id integer, tmprxn text);
\copy tmp from 'data/reaction_test_data.out.rsmi';
select * into pgreactions from (select id,reaction_from_smiles(tmprxn::cstring) rxn from tmp) as r where r is not null;
SET rdkit.move_unmmapped_reactants_to_agents=false;
select * into pgreactions_unchanged from (select id,reaction_from_smiles(tmprxn::cstring) rxn from tmp) as r where r is not null;
DROP table tmp;
SET rdkit.move_unmmapped_reactants_to_agents=true;
SELECT count(*) FROM pgreactions;
SELECT count(*) FROM pgreactions_unchanged;
SELECT SUM(reaction_numreactants(rxn)) FROM pgreactions;
SELECT SUM(reaction_numreactants(rxn)) FROM pgreactions_unchanged;
SELECT SUM(reaction_numproducts(rxn)) FROM pgreactions;
SELECT SUM(reaction_numproducts(rxn)) FROM pgreactions_unchanged;
SELECT SUM(reaction_numagents(rxn)) FROM pgreactions;
SELECT SUM(reaction_numagents(rxn)) FROM pgreactions_unchanged;
CREATE INDEX rxnidx ON pgreactions USING gist(rxn);
SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;
SELECT count(*) FROM pgreactions WHERE rxn@>'c1ccccc1>>c1ccncc1';
SELECT count(*) FROM pgreactions WHERE rxn@>'c1cnccc1>>c1ccccc1';
SELECT count(*) FROM pgreactions WHERE 'c1ccccc1>>c1ccncc1'<@rxn;
SELECT count(*) FROM pgreactions WHERE 'c1cnccc1>>c1ccccc1'<@rxn;
SELECT count(*) FROM pgreactions WHERE rxn@>'c1ccccc1>>';
SELECT count(*) FROM pgreactions WHERE rxn@>'c1cnccc1>>';
SELECT count(*) FROM pgreactions WHERE 'c1ccccc1>>'<@rxn;
SELECT count(*) FROM pgreactions WHERE 'c1cnccc1>>'<@rxn;
SELECT count(*) FROM pgreactions WHERE rxn@>'>>c1ccncc1';
SELECT count(*) FROM pgreactions WHERE rxn@>'>>c1ccccc1';
SELECT count(*) FROM pgreactions WHERE '>>c1ccncc1'<@rxn;
SELECT count(*) FROM pgreactions WHERE '>>c1ccccc1'<@rxn;
SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=off;
SELECT count(*) FROM pgreactions WHERE rxn@>'c1ccccc1>>c1ccncc1';
SELECT count(*) FROM pgreactions WHERE rxn@>'c1cnccc1>>c1ccccc1';
SELECT count(*) FROM pgreactions WHERE 'c1ccccc1>>c1ccncc1'<@rxn;
SELECT count(*) FROM pgreactions WHERE 'c1cnccc1>>c1ccccc1'<@rxn;
SELECT count(*) FROM pgreactions WHERE rxn@>'c1ccccc1>>';
SELECT count(*) FROM pgreactions WHERE rxn@>'c1cnccc1>>';
SELECT count(*) FROM pgreactions WHERE 'c1ccccc1>>'<@rxn;
SELECT count(*) FROM pgreactions WHERE 'c1cnccc1>>'<@rxn;
SELECT count(*) FROM pgreactions WHERE rxn@>'>>c1ccncc1';
SELECT count(*) FROM pgreactions WHERE rxn@>'>>c1ccccc1';
SELECT count(*) FROM pgreactions WHERE '>>c1ccncc1'<@rxn;
SELECT count(*) FROM pgreactions WHERE '>>c1ccccc1'<@rxn;
SELECT count(*) FROM pgreactions WHERE rxn?>'c1ccccc1>>c1ccncc1';
SELECT count(*) FROM pgreactions WHERE rxn?>'c1cnccc1>>c1ccccc1';
SELECT count(*) FROM pgreactions WHERE 'c1ccccc1>>c1ccncc1'?<rxn;
SELECT count(*) FROM pgreactions WHERE 'c1cnccc1>>c1ccccc1'?<rxn;
SELECT count(*) FROM pgreactions WHERE rxn?>'c1ccccc1>>';
SELECT count(*) FROM pgreactions WHERE rxn?>'c1cnccc1>>';
SELECT count(*) FROM pgreactions WHERE 'c1ccccc1>>'?<rxn;
SELECT count(*) FROM pgreactions WHERE 'c1cnccc1>>'?<rxn;
SELECT count(*) FROM pgreactions WHERE rxn?>'>>c1ccncc1';
SELECT count(*) FROM pgreactions WHERE rxn?>'>>c1ccccc1';
SELECT count(*) FROM pgreactions WHERE '>>c1ccncc1'?<rxn;
SELECT count(*) FROM pgreactions WHERE '>>c1ccccc1'?<rxn;

SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',1), reaction_difference_fp('c1ccccc1>>c1ccncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',1), reaction_difference_fp('c1ncccc1>>c1ncncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1), reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1), reaction_difference_fp('c1ncccc1>[Na+]>c1ncncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',2), reaction_difference_fp('c1ccccc1>>c1ccncc1',2));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',2), reaction_difference_fp('c1ncccc1>>c1ncncc1',2));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2), reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2), reaction_difference_fp('c1ncccc1>[Na+]>c1ncncc1',2));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',3), reaction_difference_fp('c1ccccc1>>c1ccncc1',3));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',3), reaction_difference_fp('c1ncccc1>>c1ncncc1',3));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3), reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3), reaction_difference_fp('c1ncccc1>[Na+]>c1ncncc1',3));

SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',1), reaction_structural_bfp('c1ccccc1>>c1ccncc1',1));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',1), reaction_structural_bfp('c1ncccc1>>c1ncncc1',1));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',1));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',2), reaction_structural_bfp('c1ccccc1>>c1ccncc1',2));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',2), reaction_structural_bfp('c1ncccc1>>c1ncncc1',2));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',2));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',3), reaction_structural_bfp('c1ccccc1>>c1ccncc1',3));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',3), reaction_structural_bfp('c1ncccc1>>c1ncncc1',3));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',3));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',4), reaction_structural_bfp('c1ccccc1>>c1ccncc1',4));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',4), reaction_structural_bfp('c1ncccc1>>c1ncncc1',4));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',4), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',4));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',4), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',4));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',5), reaction_structural_bfp('c1ccccc1>>c1ccncc1',5));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',5), reaction_structural_bfp('c1ncccc1>>c1ncncc1',5));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',5), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',5));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',5), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',5));

SET rdkit.agent_FP_bit_ratio=0.5;
SET rdkit.difference_FP_weight_agents=-3;
SET rdkit.difference_FP_weight_nonagents=7;

SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',1), reaction_difference_fp('c1ccccc1>>c1ccncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',1), reaction_difference_fp('c1ncccc1>>c1ncncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1), reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1), reaction_difference_fp('c1ncccc1>[Na+]>c1ncncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',2), reaction_difference_fp('c1ccccc1>>c1ccncc1',2));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',2), reaction_difference_fp('c1ncccc1>>c1ncncc1',2));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2), reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2), reaction_difference_fp('c1ncccc1>[Na+]>c1ncncc1',2));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',3), reaction_difference_fp('c1ccccc1>>c1ccncc1',3));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',3), reaction_difference_fp('c1ncccc1>>c1ncncc1',3));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3), reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3), reaction_difference_fp('c1ncccc1>[Na+]>c1ncncc1',3));

SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',1), reaction_structural_bfp('c1ccccc1>>c1ccncc1',1));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',1), reaction_structural_bfp('c1ncccc1>>c1ncncc1',1));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',1));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',2), reaction_structural_bfp('c1ccccc1>>c1ccncc1',2));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',2), reaction_structural_bfp('c1ncccc1>>c1ncncc1',2));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',2), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',2));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',3), reaction_structural_bfp('c1ccccc1>>c1ccncc1',3));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',3), reaction_structural_bfp('c1ncccc1>>c1ncncc1',3));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',3), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',3));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',4), reaction_structural_bfp('c1ccccc1>>c1ccncc1',4));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',4), reaction_structural_bfp('c1ncccc1>>c1ncncc1',4));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',4), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',4));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',4), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',4));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',5), reaction_structural_bfp('c1ccccc1>>c1ccncc1',5));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',5), reaction_structural_bfp('c1ncccc1>>c1ncncc1',5));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',5), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',5));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',5), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',5));

SET rdkit.ignore_reaction_agents=true;

SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',1), reaction_difference_fp('c1ccccc1>>c1ccncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>>c1ccncc1',1), reaction_difference_fp('c1ncccc1>>c1ncncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1), reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1));
SELECT tanimoto_sml(reaction_difference_fp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',1), reaction_difference_fp('c1ncccc1>[Na+]>c1ncncc1',1));

SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',5), reaction_structural_bfp('c1ccccc1>>c1ccncc1',5));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>>c1ccncc1',5), reaction_structural_bfp('c1ncccc1>>c1ncncc1',5));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',5), reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',5));
SELECT tanimoto_sml(reaction_structural_bfp('c1ccccc1>CC(=O)O.[Na+]>c1ccncc1',5), reaction_structural_bfp('c1ncccc1>[Na+]>c1ncncc1',5));

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=on;
DROP INDEX rxnidx;
