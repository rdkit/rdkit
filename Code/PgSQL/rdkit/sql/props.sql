SELECT mol_amw('c1ccccc1'::mol) mol_amw;
SELECT mol_logp('c1ccccc1'::mol) mol_logp;
SELECT mol_hba('c1ccccc1'::mol) mol_hba;
SELECT mol_hbd('c1ccccc1'::mol) mol_hbd;
SELECT mol_hba('c1ccncc1'::mol) mol_hba;
SELECT mol_hbd('c1ccncc1'::mol) mol_hbd;
SELECT mol_hbd('c1ccncc1O'::mol) mol_hbd;
SELECT mol_hba('c1ccncc1O'::mol) mol_hba;
SELECT mol_logp('c1ccncc1O'::mol) mol_logp;
SELECT mol_chi0n('c1ccccc1O'::mol) v;
SELECT mol_chi1n('c1ccccc1O'::mol) v;
SELECT mol_chi2n('c1ccccc1O'::mol) v;
SELECT mol_chi3n('c1ccccc1O'::mol) v;
SELECT mol_chi4n('c1ccccc1O'::mol) v;
SELECT mol_chi0v('c1ccccc1O'::mol) v;
SELECT mol_chi1v('c1ccccc1O'::mol) v;
SELECT mol_chi2v('c1ccccc1O'::mol) v;
SELECT mol_chi3v('c1ccccc1O'::mol) v;
SELECT mol_chi4v('c1ccccc1O'::mol) v;
SELECT mol_kappa1('C12CC2C3CC13'::mol) v;
SELECT mol_kappa2('CC(C)C1CCC(C)CCC1'::mol) v;
SELECT mol_kappa3('CC(C)C1CCC(C)CCC1'::mol) v;
SELECT mol_numspiroatoms('C1CCC2(C1)CC1CCC2CC1'::mol) v;
SELECT mol_numbridgeheadatoms('C1CCC2(C1)CC1CCC2CC1'::mol) v;
SELECT mol_numspiroatoms('CC1(C)CC2(C)CCC1(C)CC2'::mol) v;
SELECT mol_numbridgeheadatoms('CC1(C)CC2(C)CCC1(C)CC2'::mol) v;


-- Mol formula tests - SQL equivalents of tests in testMolDescriptors.py.
select mol_formula('[2H]C([3H])O'::mol);
                                     -- separateIsotopes = true
select mol_formula('[2H]C([3H])O'::mol, true);
                                           -- abbreviateHIsotopes = false
select mol_formula('[2H]C([3H])O'::mol, true, false);
       --
select mol_formula('[2H][13CH2]CO'::mol);
select mol_formula('[2H][13CH2]CO'::mol, true);
select mol_formula('[2H][13CH2]CO'::mol, true, false);
--
SELECT mol_numrotatablebonds('CCC'::mol) mol_numrotatablebonds;
SELECT mol_numrotatablebonds('CCCC'::mol) mol_numrotatablebonds;
SELECT mol_numrotatablebonds('c1ccccc1c1ccc(CCC)cc1'::mol) mol_numrotatablebonds;
SELECT mol_numheavyatoms('CCC'::mol) val;
SELECT mol_numatoms('CCC'::mol) val;
SELECT mol_numheteroatoms('CCC'::mol) val;
SELECT mol_numheteroatoms('CCO'::mol) val;
SELECT mol_tpsa('CCC'::mol) val;
SELECT mol_tpsa('CCO'::mol) val;
SELECT mol_numrings('CCC'::mol) val;
SELECT mol_numrings('C1CC1'::mol) val;
SELECT mol_murckoscaffold('c1ccccc1CCC'::mol) val;
SELECT mol_murckoscaffold('CSC(C)=O'::mol) is NULL;

SELECT mol_to_svg('CCO'::mol) svg;
SELECT mol_to_svg('CCO'::mol,'legend') svg;
SELECT mol_to_svg('CCO'::mol,'legend',250,200,
  '{"atomLabels":{"1":"foo"},"legendColour":[0.5,0.5,0.5]}') svg;

-- GitHub Issue 2174 - mol_to_svg() should not change input mol.
/**
  Check that mol_to_svg() does not change the_mol.
  In previous versions, the atom+bond count would
  change from '10 11' in "before_mol" to '11 12'
  in "after_mol", due to mol_to_svg()'s call to
  MolDraw2DUtils::prepareMolForDrawing().
**/
with t as (
  select 'C[C@H]1CC[C@H]2CCCCC12'::mol as the_mol
)
select
    substring(mol_to_ctab(the_mol)::text, 1, 65) as before_mol
  , substring(mol_to_svg(the_mol)::text, 1, 20) as the_svg
  , substring(mol_to_ctab(the_mol)::text, 1, 65) as after_mol
  from t;
