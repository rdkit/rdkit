//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

package org.RDKit;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.List;

import org.junit.*;

public class QueryAtomTest extends GraphMolTest {

    @Test
    public void testChemdrawBlock() {
        // smiles c1[nH]ccc1
        // Chemdraw will convert explicit hydrogens in aromatic compounds to query hydrogen count
        var block = String.join("\n",
                "ChemDraw05012418152D",
                "",
                "0 0 0 0 0 0 V3000",
                "M V30 BEGIN CTAB",
                "M V30 COUNTS 5 5 0 0 0",
                "M V30 BEGIN ATOM",
                "M V30 1 C - 0.652854 - 0.050688 0.000000 0",
                "M V30 2 C - 0.240354 0.663783 0.000000 0",
                "M V30 3 C 0.566618 0.492256 0.000000 0",
                "M V30 4 C 0.652854 - 0.328225 0.000000 0",
                "M V30 5 N - 0.100821 - 0.663783 0.000000 0 HCOUNT = 1",
                "M V30 END ATOM",
                "M V30 BEGIN BOND",
                "M V30 1 4 1 2",
                "M V30 2 4 2 3",
                "M V30 3 4 3 4",
                "M V30 4 4 4 5",
                "M V30 5 4 1 5",
                "M V30 END BOND",
                "M V30 END CTAB",
                "M END",
                "");
        RWMol mol = RWMol.MolFromMolBlock(block, false);
        assertNotNull(mol);
        assertEquals(5, mol.getNumAtoms());
        var queryAtom = mol.getAtomWithIdx(4);
        assertTrue(queryAtom.hasQuery());
        var smarts = queryAtom.AtomGetSmarts();
        var expected = "[#7&h{1-}]";
        assertEquals(expected, smarts);
        var description = queryAtom.describeQuery();
        expected = "AtomAnd\n  AtomAtomicNum 7 = val\n  less_AtomImplicitHCount 1 <= \n";
        assertEquals(expected, description);

        var molSmiles = mol.MolToSmiles();

        var atom = new Atom(7);
        atom.setNumExplicitHs(1);
        mol.replaceAtom(4, atom);
        mol.clearComputedProps(true);
        mol.sanitizeMol();

        atom = mol.getAtomWithIdx(4);
        assertEquals(7, atom.getAtomicNum());
        assertEquals(1, atom.getNumExplicitHs());
        molSmiles = mol.MolToSmiles();
        assertEquals("c1cc[nH]c1", molSmiles);
    }

    @Test
    public void testReplaceQuerySimple() {
        var query = RWMol.MolFromSmiles("c1[nH]ccc1");
        var nitrogen = query.getAtomWithIdx(2);
        var imp = nitrogen.getNumImplicitHs();
        assertEquals(1, imp);
        var queryAtom = RDKFuncs.replaceAtomWithQueryAtom(query, nitrogen);
        queryAtom.ExpandQuery(RDKFuncs.ExplicitDegreeEqualsQueryAtom(3));
        var mol1 = RWMol.MolFromSmiles("Cc1[nH]ccc1");
        assertTrue(mol1.hasSubstructMatch(query));
        var mol2 = RWMol.MolFromSmiles("c1[nH]ccc1");
        assertFalse(mol2.hasSubstructMatch(query));
    }

    @Test
    public void testExpandQuerySimple() {
        var query = RWMol.MolFromSmarts("c1nccc1");
        var mol1 = RWMol.MolFromSmiles("c1[nH]ccc1");
        assertTrue(mol1.hasSubstructMatch(query));
        var atom = query.getAtomWithIdx(0);
        assertTrue(atom.hasQuery());
        var queryAtom = RDKFuncs.ExplicitDegreeEqualsQueryAtom(3);
        atom.ExpandQuery(queryAtom, CompositeQueryType.COMPOSITE_AND);
        assertFalse(mol1.hasSubstructMatch(query));
        var mol2 = RWMol.MolFromSmiles("Cc1[nH]ccc1");
        assertTrue(mol2.hasSubstructMatch(query));
    }

    // adapted from Python tests

    private int[] GetAtomIdxsMatchingQuery(ROMol mol, QueryAtom qa) {
        List<Integer> indices = new ArrayList<Integer>();
        for (int i =0; i< mol.getNumAtoms(); i++) {
            var atom = mol.getAtomWithIdx(i);
            if (qa.MatchAtom(atom)) {
                indices.add(i);
            }
        }
        int[] arr = new int[indices.size()];
        for (int i=0; i<arr.length; i++) {
            arr[i] = indices.get(i);
        }
        return arr;
    }

    @Test
    public void testQueryAtoms() {
        var m = RWMol.MolFromSmiles("c1nc(C)n(CC)c1");
        var qa = RDKFuncs.ExplicitDegreeEqualsQueryAtom(3);
        var matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {2, 4}, matches);

        qa.ExpandQuery(RDKFuncs.AtomNumEqualsQueryAtom(6, true));
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {4}, matches);

        qa = RDKFuncs.ExplicitDegreeEqualsQueryAtom(3);
        qa.ExpandQuery(RDKFuncs.AtomNumEqualsQueryAtom(6, true),
                CompositeQueryType.COMPOSITE_OR);
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {1, 2, 4}, matches);

        qa = RDKFuncs.ExplicitDegreeEqualsQueryAtom(3);
        qa.ExpandQuery(RDKFuncs.AtomNumEqualsQueryAtom(6, true),
                CompositeQueryType.COMPOSITE_XOR);
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {1, 2}, matches);

        qa = RDKFuncs.ExplicitDegreeGreaterQueryAtom(2);
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {2, 4}, matches);

        qa = RDKFuncs.ExplicitDegreeLessQueryAtom(2);
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {3, 6}, matches);

        m = RWMol.MolFromSmiles("N[CH][CH]");
        qa = RDKFuncs.NumRadicalElectronsGreaterQueryAtom(0);
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {1, 2}, matches);
        qa = RDKFuncs.NumRadicalElectronsGreaterQueryAtom(1);
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {2}, matches);

        m = RWMol.MolFromSmiles("F[C@H](Cl)C");
        qa = RDKFuncs.HasChiralTagQueryAtom();
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {1}, matches);
        qa = RDKFuncs.MissingChiralTagQueryAtom();
        assertArrayEquals(new int[] {1}, matches);

        m = RWMol.MolFromSmiles("F[CH](Cl)C");
        qa = RDKFuncs.HasChiralTagQueryAtom();
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {}, matches);
        qa = RDKFuncs.MissingChiralTagQueryAtom();
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {1}, matches);

        m = RWMol.MolFromSmiles("CNCON");
        qa = RDKFuncs.NumHeteroatomNeighborsEqualsQueryAtom(2);
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {2}, matches);
        qa = RDKFuncs.NumHeteroatomNeighborsGreaterQueryAtom(0);
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {0, 2, 3, 4}, matches);

        m = RWMol.MolFromSmiles("CC12CCN(CC1)C2");
        qa = RDKFuncs.IsBridgeheadQueryAtom();
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {1, 4}, matches);

        m = RWMol.MolFromSmiles("OCCOC");
        qa = RDKFuncs.NonHydrogenDegreeEqualsQueryAtom(2);
        matches = GetAtomIdxsMatchingQuery(m, qa);
        assertArrayEquals(new int[] {1, 2, 3}, matches);
    }

    private int[] GetBondIdxsMatchingQuery(ROMol mol, QueryBond qa) {
        List<Integer> indices = new ArrayList<Integer>();
        for (int i =0; i< mol.getNumBonds(); i++) {
            var atom = mol.getBondWithIdx(i);
            if (qa.MatchBond(atom)) {
                indices.add(i);
            }
        }
        int[] arr = new int[indices.size()];
        for (int i=0; i<arr.length; i++) {
            arr[i] = indices.get(i);
        }
        return arr;
    }

    @Test
    public void TestBondPropQueries() {
        var m = RWMol.MolFromSmiles("CCCCCCCCCCCCCC");
        var bonds = m.getBonds();
        bonds.get(0).setProp("hah", "hah");
        bonds.get(1).setIntProp("bar", 1);
        bonds.get(2).setIntProp("bar", 2);
        bonds.get(3).setBoolProp("baz", true);
        bonds.get(4).setBoolProp("baz", false);
        bonds.get(5).setProp("boo", "hoo");
        bonds.get(6).setProp("boo", "-urns");
        bonds.get(7).setDoubleProp("boot", 1.0);
        bonds.get(8).setDoubleProp("boot", 4.0);
        bonds.get(9).setDoubleProp("number", 4.0);
        bonds.get(10).setIntProp("number", 4);

        var qb = RDKFuncs.HasIntPropWithValueQueryBond("bar", 1);
        var matches = GetBondIdxsMatchingQuery(m, qb);
        assertArrayEquals(new int[] {1}, matches);

        qb = RDKFuncs.HasIntPropWithValueQueryBond("bar", 2);
        matches = GetBondIdxsMatchingQuery(m, qb);
        assertArrayEquals(new int[] {2}, matches);

        qb = RDKFuncs.HasBoolPropWithValueQueryBond("baz", true);
        matches = GetBondIdxsMatchingQuery(m, qb);
        assertArrayEquals(new int[] {3}, matches);

        qb = RDKFuncs.HasBoolPropWithValueQueryBond("baz", false);
        matches = GetBondIdxsMatchingQuery(m, qb);
        assertArrayEquals(new int[] {4}, matches);

        qb = RDKFuncs.HasStringPropWithValueQueryBond("boo", "hoo");
        matches = GetBondIdxsMatchingQuery(m, qb);
        assertArrayEquals(new int[] {5}, matches);

        qb = RDKFuncs.HasStringPropWithValueQueryBond("boo", "-urns");
        matches = GetBondIdxsMatchingQuery(m, qb);
        assertArrayEquals(new int[] {6}, matches);

        qb = RDKFuncs.HasDoublePropWithValueQueryBond("boot", 1.0);
        matches = GetBondIdxsMatchingQuery(m, qb);
        assertArrayEquals(new int[] {7}, matches);

        qb = RDKFuncs.HasDoublePropWithValueQueryBond("boot", 4.0);
        matches = GetBondIdxsMatchingQuery(m, qb);
        assertArrayEquals(new int[] {8}, matches);

        qb = RDKFuncs.HasDoublePropWithValueQueryBond("boot", 1.0, false, 3.0);
        matches = GetBondIdxsMatchingQuery(m, qb);
        assertArrayEquals(new int[] {7, 8}, matches);

        qb = RDKFuncs.HasIntPropWithValueQueryBond("number", 4);
        matches = GetBondIdxsMatchingQuery(m, qb);
        assertArrayEquals(new int[] {10}, matches);
    }

    public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("org.RDKit.RascalMCESTest");
    }
}
