
import doctest
import unittest

import os
import csv
from random import Random
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateHeterocycles import (
    GetHeterocycleReactionSmarts,
    GetHeterocycleReactions,
    EnumerateHeterocycles,
)

def has_radical(mol):
    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons():
            return True
    return False

def has_aromatic(mol):
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            return True
    return False

def check_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    assert not has_radical(mol)
    assert len(Chem.rdmolops.GetMolFrags(mol)) == 1, smiles
    return mol

def assert_valid_change(orig_can_smi, smiles):
    orig_mol = Chem.MolFromSmiles(orig_can_smi)
    assert orig_can_smi != smiles, "%s != %s" % (orig_can_smi, smiles)
    msg = "%s produced %s" % (orig_can_smi, smiles)
    try:
        mol = check_smiles(smiles)
    except Exception:
        print(msg)
        raise

    assert has_aromatic(mol), 'lost aromaticity when ' + msg

    orig_charge = Chem.rdmolops.GetFormalCharge(orig_mol)
    new_charge = Chem.rdmolops.GetFormalCharge(mol)

    if orig_charge != new_charge:
        assert orig_charge - 1 == new_charge, msg

def get_unique_products(rxn, mol):
    uniq_smiles = set()
    for newmol in rxn.RunReactants((mol,)):
        isosmiles = Chem.MolToSmiles(newmol[0], isomericSmiles=True)
        if isosmiles in uniq_smiles:
            continue
        uniq_smiles.add(isosmiles)
        newmol = Chem.MolFromSmiles(isosmiles)
        assert newmol is not None, '%s produced %s when applying' % (Chem.MolToSmiles(mol, isomericSmiles=True),
                                                                     isosmiles)
        yield Chem.MolToSmiles(newmol, isomericSmiles=True)


class TestCase(unittest.TestCase):
    def test_smarts_match_examples(self):
        for row in GetHeterocycleReactionSmarts():
            smarts = row.SMARTS
            if not smarts:
                continue

            substruct = Chem.MolFromSmarts(smarts)
            for smiles in row.EXAMPLE.split(','):
                assert smiles
                mol = Chem.MolFromSmiles(smiles)
                assert mol.HasSubstructMatch(substruct), "%s not in %s" % (smarts, smiles)

            for smiles in row.NEGATIVE_EXAMPLE.split(','):
                if not smiles:
                    continue

                mol = Chem.MolFromSmiles(smiles)
                assert not mol.HasSubstructMatch(substruct), "%s should not be in %s" % (smarts, smiles)

    def test_at_least_something_matches_every_negative_example(self):
        negative_examples = []
        substructs = []
        for row in GetHeterocycleReactionSmarts():
            for smiles in row.NEGATIVE_EXAMPLE.split(','):
                if not smiles:
                    continue
                mol = Chem.MolFromSmiles(smiles)
                assert mol is not None
                negative_examples.append(mol)

            smarts = row.SMARTS
            if not smarts:
                continue

            substruct = Chem.MolFromSmarts(smarts)
            assert substruct is not None
            substructs.append(substruct)

        for mol in negative_examples:
            something_hit = False
            for substruct in substructs:
                if mol.HasSubstructMatch(substruct):
                    something_hit = True
                    break
            assert something_hit, ('nothing matched %s' % Chem.MolToSmiles(mol, isomericSmiles=True))

    def test_reactions_modify_examples(self):
        for row in GetHeterocycleReactionSmarts():
            smarts = row.SMARTS
            if not smarts:
                continue

            for product in row.CONVERT_TO.split(','):
                reaction = smarts + '>>' + product
                rxn = AllChem.ReactionFromSmarts(reaction)

                for smiles in row.EXAMPLE.split(','):
                    orig_can_smi = Chem.CanonSmiles(smiles)
                    assert smiles
                    mol = Chem.MolFromSmiles(smiles)
                    for newmol in rxn.RunReactants((mol,)):
                        newmol = newmol[0]
                        isosmi = Chem.MolToSmiles(newmol, isomericSmiles=True)
                        assert_valid_change(orig_can_smi, isosmi)

    def test_apply_every_rule_to_every_fragment(self):
        fieldnames = ['SMILES', 'MUTATED', 'REACTION', 'DESCRIPTION']
        writer = csv.DictWriter(open('hetero_atom_mutations.csv', 'w'), fieldnames)
        writer.writeheader()

        notchanged = csv.DictWriter(open('not_changed.csv', 'w'), ['SMILES', 'TITLE'])
        notchanged.writeheader()

        fragment_library = os.path.join(os.path.dirname(__file__), 'test_data', 'fragments.csv')
        frag_reader = csv.DictReader(open(fragment_library))
        for row in frag_reader:
            smiles = row['SMILES']
            rdkit_mol = Chem.MolFromSmiles(smiles)
            orig_can_smi = Chem.MolToSmiles(rdkit_mol, isomericSmiles=True)

            changed = False

            for src, rxn in zip(GetHeterocycleReactionSmarts(), GetHeterocycleReactions()):
                for smiles in get_unique_products(rxn, rdkit_mol):
                    assert_valid_change(orig_can_smi, smiles)
                    row = {
                        'SMILES'  : orig_can_smi,
                        'MUTATED' : smiles,
                        'REACTION' : src.SMARTS + '>>' + src.CONVERT_TO,
                        'DESCRIPTION' : src.DESCRIPTION,
                        }
                    writer.writerow(row)
                    changed = True

            # record aromatic fragments that no rule changes (possible problems?)
            if not changed and has_aromatic(rdkit_mol):
                row = {'SMILES' : orig_can_smi,
                       'TITLE'  : orig_can_smi}
                notchanged.writerow(row)

# There are only 2 possible 6 member ring 3 heteroatom systems, 3 carbons in each (operating under the assumption we never seen 4+!)
# c1ccnnn1 => c1ccnnn1, c1cnnnc1, c1nnncc1
# c1cnncn1 => c1cnncn1, c1nncnc1, c1ncnnc1

# Or we can enumerate all carbons for every 2 heteroatom system, 10 total (safer as it if 4 heteroatoms are passed in, it won't add more)
# c1cccnn1 => c1cccnn1, c1ccnnc1, c1cnncc1, c1nnccc1
# c1ccncn1 => c1ccncn1, c1cncnc1, c1ncncc1, c1ncccn1
# c1cnccn1 => c1cnccn1, c1nccnc1 (other 2 are symmetry)

# $(a1cccnn1),$(a1ccnnc1),$(a1cnncc1),$(a1nnccc1),$(a1ccncn1),$(a1cncnc1),$(a1ncncc1),$(a1ncccn1),$(a1cnccn1),$(a1nccnc1) =>
# only one row needs to allow carbons
# [c;h1;D2;r6;$(a1ccc[c,n,o][c,n,o]1),$(a1cc[c,n,o][c,n,o]c1),$(a1c[c,n,o][c,n,o]cc1),$(a1[c,n,o][c,n,o]ccc1),$(a1cc[n,o]c[n,o]1),$(a1c[n,o]c[n,o]c1),$(a1[n,o]c[n,o]cc1),$(a1[n,o]ccc[n,o]1),$(a1c[n,o]cc[n,o]1),$(a1[n,o]cc[n,o]c1):1]

    def get_six_member_ring_carbon_to_nitrogen_reaction(self):
        expected_description = 'aromatic carbon in 6 membered ring'
        rxns = [r for r in GetHeterocycleReactionSmarts() if r.DESCRIPTION.startswith(expected_description)]
        assert len(rxns) == 1, "expecting only one of these rules for now"
        return rxns[0]

    def test_six_member_ring_carbon_to_nitrogen_should_hit_all_other_carbons(self):
        rxn = self.get_six_member_ring_carbon_to_nitrogen_reaction()

        test_inputs = [
            #SMILES, num expected matches
            ('c1ccccc1', 6),
            ('n1ccccc1', 5),
            ('n1ncccc1', 4),
            ('n1cnccc1', 4),
            ('n1ccncc1', 4),
            ('c1c(*)[nH]c(=O)[nH]c1=O', 1),
            ('*c1c[nH]c(=O)[nH]c1=O', 1),
            ]

        smarts_mol = Chem.MolFromSmarts(rxn.SMARTS)

        for smiles, expected_num_matches in test_inputs:
            rdkit_mol = Chem.MolFromSmiles(smiles)
            num_rdkit_matches = len(rdkit_mol.GetSubstructMatches(smarts_mol))

            assert num_rdkit_matches == expected_num_matches, "expected %u matches in %s, got %u" % (expected_num_matches, smiles, num_rdkit_matches)


    def test_fuzz_atom_mutations(self):
        fragment_library = os.path.join(os.path.dirname(__file__), 'test_data', 'fragments.csv')
        base, ext = os.path.splitext(os.path.basename(fragment_library))

        rand = Random(0xDEADBEEF)
        uniq_fragments = set()
        fragments = []

        frag_reader = csv.DictReader(open(fragment_library))
        for row in frag_reader:
            smiles = row['SMILES']
            rdkit_mol = Chem.MolFromSmiles(smiles)
            if not has_aromatic(rdkit_mol):
                continue
            orig_can_smi = Chem.MolToSmiles(rdkit_mol, isomericSmiles=True)
            assert orig_can_smi not in uniq_fragments
            uniq_fragments.add(orig_can_smi)
            fragments.append(orig_can_smi)

        print(len(fragments), "fragments with aromaticity")

        fieldnames = ['SMILES', 'MUTATED', 'REACTION', 'DESCRIPTION']
        writer = csv.DictWriter(open(base + 'hetero_atom_mutations_fuzzing.csv', 'w'), fieldnames)
        writer.writeheader()
        notchanged = csv.DictWriter(open(base + 'not_changed_during_fuzzing.csv', 'w'), ['SMILES', 'TITLE'])
        notchanged.writeheader()
        uniq_notchanged = set()

        num_trials = 1000
        # to test the full range of possible fragments
        #num_trials = 1000000

        total_generated = 0
        for i in range(num_trials):
            if i and i % 1000 == 0:
                print(i)
            if not fragments:
                print("Converged! No more fragments left!")
                break
            idx = rand.randint(0, len(fragments) - 1)
            orig_can_smi = fragments.pop(idx)
            rdkit_mol = Chem.MolFromSmiles(orig_can_smi)

            changed = False

            for src, rxn in zip(GetHeterocycleReactionSmarts(), GetHeterocycleReactions()):
                for smiles in get_unique_products(rxn, rdkit_mol):
                    total_generated += 1
                    changed = True
                    assert_valid_change(orig_can_smi, smiles)
                    if smiles in uniq_fragments:
                        continue
                    uniq_fragments.add(smiles)
                    fragments.append(smiles)

                    row = {
                        'SMILES'  : orig_can_smi,
                        'MUTATED' : smiles,
                        'REACTION' : src.SMARTS + '>>' + src.CONVERT_TO,
                        'DESCRIPTION' : src.DESCRIPTION,
                        }
                    writer.writerow(row)


            # record aromatic fragments that no rule changes (possible problems?)
            if not changed and orig_can_smi not in uniq_notchanged:
                uniq_notchanged.add(orig_can_smi)
                row = {'SMILES' : orig_can_smi,
                       'TITLE'  : orig_can_smi}
                notchanged.writerow(row)
        print(total_generated, "generated of which", len(uniq_fragments), "are unique fragments generated after", num_trials, "trials")

    def test_reaction_on_mol_with_attachments(self):
        src = self.get_six_member_ring_carbon_to_nitrogen_reaction()
        reaction = src.SMARTS + '>>' + src.CONVERT_TO
        rxn = AllChem.ReactionFromSmarts(reaction)


        rdkit_mol = Chem.MolFromSmiles('c1([*:2])ccccc1([*:3])')
        products = set(get_unique_products(rxn, rdkit_mol))
        expected = set(['c1cnc([*:2])c([*:3])c1',
                        'c1cc([*:3])c([*:2])cn1',
                        'c1cnc([*:3])c([*:2])c1',
                        'c1cc([*:2])c([*:3])cn1'])
        assert products == expected, '%r != %r' % (products, expected)

        # two less of these because of symmetry
        rdkit_mol = Chem.MolFromSmiles('c1([*:2])ccc([*:3])cc1')
        assert set(get_unique_products(rxn, rdkit_mol)) == set(['c1cc([*:2])ncc1[*:3]',
                                                                'c1cc([*:3])ncc1[*:2]'])

        # test what happens when nothing matches
        rdkit_mol = Chem.MolFromSmiles('c1([*:2])cc[nH]c1([*:3])')
        assert set(get_unique_products(rxn, rdkit_mol)) == set()


    def test_enumerate_heterocycles(self):
        rand = Random(0xBAD1DEA)

        fragment_library = os.path.join(os.path.dirname(__file__), 'test_data', 'fragments.csv')
        frag_reader = csv.DictReader(open(fragment_library))
        for row in frag_reader:
            smiles = row['SMILES']
            mol = Chem.MolFromSmiles(smiles)
            reference_set = set(Chem.MolToSmiles(m) for m in EnumerateHeterocycles(mol))

            if not has_aromatic(Chem.MolFromSmiles(smiles)):
                assert len(reference_set) == 0
                continue

            assert len(reference_set) > 0

            for nbr_smiles in rand.sample(list(reference_set), 4):
                nbr = Chem.MolFromSmiles(nbr_smiles)
                comparison_set = set(Chem.MolToSmiles(m) for m in EnumerateHeterocycles(nbr))
                assert reference_set == comparison_set


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
