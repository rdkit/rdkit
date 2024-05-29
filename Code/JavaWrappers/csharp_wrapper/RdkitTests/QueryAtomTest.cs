//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.PortableExecutable;
using GraphMolWrap;
using Xunit;

namespace RdkitTests
{
    public class QueryAtomTest
    {
        [Fact]
        public void TestChemdrawBlock()
        {
            // smiles c1[nH]ccc1
            // Chemdraw will convert explicit hydrogens in aromatic compounds to query hydrogen count
            var block = @"
  ChemDraw05012418152D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 5 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.652854 -0.050688 0.000000 0
M  V30 2 C -0.240354 0.663783 0.000000 0
M  V30 3 C 0.566618 0.492256 0.000000 0
M  V30 4 C 0.652854 -0.328225 0.000000 0
M  V30 5 N -0.100821 -0.663783 0.000000 0 HCOUNT=1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 4 1 2
M  V30 2 4 2 3
M  V30 3 4 3 4
M  V30 4 4 4 5
M  V30 5 4 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
";
            var mol = RWMol.MolFromMolBlock(block, false);
            Assert.NotNull(mol);
            Assert.Equal(5U, mol.getNumAtoms());
            var queryAtom = mol.getAtomWithIdx(4);
            Assert.True(queryAtom.hasQuery());
            var smarts = queryAtom.AtomGetSmarts();
            var expected = "[#7&h{1-}]";
            Assert.Equal(expected, smarts);
            var description = queryAtom.describeQuery();
            expected = "AtomAnd\n  AtomAtomicNum 7 = val\n  less_AtomImplicitHCount 1 <= \n";
            Assert.Equal(expected, description);

            var molSmiles = mol.MolToSmiles();
            var molSmarts = RDKFuncs.MolToSmarts(mol);

            var atom = new Atom(7);
            atom.setNumExplicitHs(1);
            mol.replaceAtom(4, atom);
            mol.clearComputedProps(true);
            mol.sanitizeMol();

            atom = mol.getAtomWithIdx(4);
            Assert.Equal(7, atom.getAtomicNum());
            Assert.Equal(1U, atom.getNumExplicitHs());
            molSmiles = mol.MolToSmiles();
            Assert.Equal("c1cc[nH]c1", molSmiles);
        }

        [Fact]
        public void TestReplaceQuerySimple()
        {
            var query = RWMol.MolFromSmiles("c1[nH]ccc1");
            var nitrogen = query.getAtomWithIdx(2);
            var imp = nitrogen.getNumImplicitHs();
            Assert.Equal(1U, imp);
            var queryAtom = RDKFuncs.replaceAtomWithQueryAtom(query, nitrogen);
            queryAtom.ExpandQuery(RDKFuncs.ExplicitDegreeEqualsQueryAtom(3));
            var mol1 = RWMol.MolFromSmiles("Cc1[nH]ccc1");
            Assert.True(mol1.hasSubstructMatch(query));
            var mol2 = RWMol.MolFromSmiles("c1[nH]ccc1");
            Assert.False(mol2.hasSubstructMatch(query));
        }

        [Fact]
        public void TestExpandQuerySimple()
        {
            var query = RWMol.MolFromSmarts("c1nccc1");
            var mol1 = RWMol.MolFromSmiles("c1[nH]ccc1");
            Assert.True(mol1.hasSubstructMatch(query));
            var atom = query.getAtomWithIdx(0);
            Assert.True(atom.hasQuery());
            var queryAtom = RDKFuncs.ExplicitDegreeEqualsQueryAtom(3);
            atom.ExpandQuery(queryAtom, CompositeQueryType.COMPOSITE_AND);
            Assert.False(mol1.hasSubstructMatch(query));
            var mol2 = RWMol.MolFromSmiles("Cc1[nH]ccc1");
            Assert.True(mol2.hasSubstructMatch(query));
        }


        // adapted from Python tests
        private IList<uint> GetAtomIdxsMatchingQuery(ROMol mol, QueryAtom qa)
        {
            return mol.getAtoms().Where(qa.MatchAtom).Select(a => a.getIdx()).ToList();
        }

        [Fact]
        public void TestQueryAtoms()
        {
            var m = RWMol.MolFromSmiles("c1nc(C)n(CC)c1");
            var qa = RDKFuncs.ExplicitDegreeEqualsQueryAtom(3);
            var matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 2U, 4U }, matches);

            qa.ExpandQuery(RDKFuncs.AtomNumEqualsQueryAtom(6, true));
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 4U }, matches);

            qa = RDKFuncs.ExplicitDegreeEqualsQueryAtom(3);
            qa.ExpandQuery(RDKFuncs.AtomNumEqualsQueryAtom(6, true),
                CompositeQueryType.COMPOSITE_OR);
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 1U, 2U, 4U }, matches);

            qa = RDKFuncs.ExplicitDegreeEqualsQueryAtom(3);
            qa.ExpandQuery(RDKFuncs.AtomNumEqualsQueryAtom(6, true),
                CompositeQueryType.COMPOSITE_XOR);
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 1U, 2U }, matches);

            qa = RDKFuncs.ExplicitDegreeGreaterQueryAtom(2);
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 2U, 4U }, matches);

            qa = RDKFuncs.ExplicitDegreeLessQueryAtom(2);
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 3U, 6U }, matches);

            m = RWMol.MolFromSmiles("N[CH][CH]");
            qa = RDKFuncs.NumRadicalElectronsGreaterQueryAtom(0);
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 1U, 2U }, matches);
            qa = RDKFuncs.NumRadicalElectronsGreaterQueryAtom(1);
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 2U }, matches);

            m = RWMol.MolFromSmiles("F[C@H](Cl)C");
            qa = RDKFuncs.HasChiralTagQueryAtom();
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 1U }, matches);
            qa = RDKFuncs.MissingChiralTagQueryAtom();
            Assert.Equal(new List<uint> { 1U }, matches);

            m = RWMol.MolFromSmiles("F[CH](Cl)C");
            qa = RDKFuncs.HasChiralTagQueryAtom();
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { }, matches);
            qa = RDKFuncs.MissingChiralTagQueryAtom();
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 1U }, matches);

            m = RWMol.MolFromSmiles("CNCON");
            qa = RDKFuncs.NumHeteroatomNeighborsEqualsQueryAtom(2);
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 2U }, matches);
            qa = RDKFuncs.NumHeteroatomNeighborsGreaterQueryAtom(0);
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 0U, 2U, 3U, 4U }, matches);

            m = RWMol.MolFromSmiles("CC12CCN(CC1)C2");
            qa = RDKFuncs.IsBridgeheadQueryAtom();
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 1U, 4U }, matches);

            m = RWMol.MolFromSmiles("OCCOC");
            qa = RDKFuncs.NonHydrogenDegreeEqualsQueryAtom(2);
            matches = GetAtomIdxsMatchingQuery(m, qa);
            Assert.Equal(new List<uint> { 1U, 2U, 3U }, matches);
        }
        private IList<uint> GetBondIdxsMatchingQuery(ROMol mol, QueryBond qb)
        {
            return mol.getBonds().Where(qb.MatchBond).Select(a => a.getIdx()).ToList();
        }


        [Fact]
        public void TestBondPropQueries()
        {
            var m = RWMol.MolFromSmiles("CCCCCCCCCCCCCC");
            var bonds = m.getBonds();
            bonds[0].setProp("hah", "hah");
            bonds[1].setIntProp("bar", 1);
            bonds[2].setIntProp("bar", 2);
            bonds[3].setBoolProp("baz", true);
            bonds[4].setBoolProp("baz", false);
            bonds[5].setProp("boo", "hoo");
            bonds[6].setProp("boo", "-urns");
            bonds[7].setDoubleProp("boot", 1.0);
            bonds[8].setDoubleProp("boot", 4.0);
            bonds[9].setDoubleProp("number", 4.0);
            bonds[10].setIntProp("number", 4);

            var qb = RDKFuncs.HasIntPropWithValueQueryBond("bar",1);
            var matches = GetBondIdxsMatchingQuery(m, qb);
            Assert.Equal(new List<uint> { 1U }, matches);

            qb = RDKFuncs.HasIntPropWithValueQueryBond("bar",2);
            matches = GetBondIdxsMatchingQuery(m, qb);
            Assert.Equal(new List<uint> { 2U }, matches);

            qb = RDKFuncs.HasBoolPropWithValueQueryBond("baz",true);
            matches = GetBondIdxsMatchingQuery(m, qb);
            Assert.Equal(new List<uint> { 3U }, matches);

            qb = RDKFuncs.HasBoolPropWithValueQueryBond("baz",false);
            matches = GetBondIdxsMatchingQuery(m, qb);
            Assert.Equal(new List<uint> { 4U }, matches);

            qb = RDKFuncs.HasStringPropWithValueQueryBond("boo","hoo");
            matches = GetBondIdxsMatchingQuery(m, qb);
            Assert.Equal(new List<uint> { 5U }, matches);

            qb = RDKFuncs.HasStringPropWithValueQueryBond("boo","-urns");
            matches = GetBondIdxsMatchingQuery(m, qb);
            Assert.Equal(new List<uint> { 6U }, matches);

            qb = RDKFuncs.HasDoublePropWithValueQueryBond("boot",1.0);
            matches = GetBondIdxsMatchingQuery(m, qb);
            Assert.Equal(new List<uint> { 7U }, matches);

            qb = RDKFuncs.HasDoublePropWithValueQueryBond("boot",4.0);
            matches = GetBondIdxsMatchingQuery(m, qb);
            Assert.Equal(new List<uint> { 8U }, matches);

            qb = RDKFuncs.HasDoublePropWithValueQueryBond("boot",1.0, false, 3.0);
            matches = GetBondIdxsMatchingQuery(m, qb);
            Assert.Equal(new List<uint> { 7U, 8U }, matches);

            qb = RDKFuncs.HasIntPropWithValueQueryBond("number",4);
            matches = GetBondIdxsMatchingQuery(m, qb);
            Assert.Equal(new List<uint> { 10U }, matches);
        }
    }
}