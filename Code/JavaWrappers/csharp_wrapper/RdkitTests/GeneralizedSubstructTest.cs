using System;
using GraphMolWrap;
using Xunit;

namespace RdkitTests
{
    public class GeneralizedSubstructTest
    {
        private bool FingerprintsMatch(ExtendedQueryMol queryMol, RWMol target)
        {
            var queryFingerprint = queryMol.patternFingerprintQuery();
            var targetFingerprint = RDKFuncs.patternFingerprintTargetMol(target);
            Assert.True(queryFingerprint.getNumOnBits() > 0);
            Assert.True(targetFingerprint.getNumOnBits() > 0);
            var match = RDKFuncs.AllProbeBitsMatch(queryFingerprint, targetFingerprint);
            return match;
        }

        [Fact]
        public void TestControlSteps()
        {
            var queryMol = RWMol.MolFromSmiles("COCC1OC(N)=N1 |LN:1:1.3|");
            var xqm1 = RDKFuncs.createExtendedQueryMol(queryMol);
            var xqm2 = RDKFuncs.createExtendedQueryMol(queryMol, false);
            var xqm3 = RDKFuncs.createExtendedQueryMol(queryMol, true, false);
            var xqm4 = RDKFuncs.createExtendedQueryMol(queryMol, false, false);

            var mol1 = RWMol.MolFromSmiles("COCC1OC(N)=N1");
            Assert.Equal(1, xqm1.getSubstructMatches(mol1).Count);
            Assert.Equal(1, xqm2.getSubstructMatches(mol1).Count);
            Assert.Equal(1, xqm3.getSubstructMatches(mol1).Count);
            Assert.Equal(1, xqm4.getSubstructMatches(mol1).Count);
            Assert.True(FingerprintsMatch(xqm1, mol1));
            Assert.True(FingerprintsMatch(xqm2, mol1));
            Assert.True(FingerprintsMatch(xqm3, mol1));
            Assert.True(FingerprintsMatch(xqm4, mol1));

            var mol2 = RWMol.MolFromSmiles("COCC1OC(=N)N1");
            Assert.Equal(1, xqm1.getSubstructMatches(mol2).Count);
            Assert.Equal(1, xqm2.getSubstructMatches(mol2).Count);
            Assert.Equal(0, xqm3.getSubstructMatches(mol2).Count);
            Assert.Equal(0, xqm4.getSubstructMatches(mol2).Count);
            Assert.True(FingerprintsMatch(xqm1, mol2));
            Assert.True(FingerprintsMatch(xqm2, mol2));
            Assert.False(FingerprintsMatch(xqm3, mol2));
            Assert.False(FingerprintsMatch(xqm4, mol2));

            var mol3 = RWMol.MolFromSmiles("COOCC1OC(N)=N1");
            Assert.Equal(1, xqm1.getSubstructMatches(mol3).Count);
            Assert.Equal(0, xqm2.getSubstructMatches(mol3).Count);
            Assert.Equal(1, xqm3.getSubstructMatches(mol3).Count);
            Assert.Equal(0, xqm4.getSubstructMatches(mol3).Count);
            Assert.True(FingerprintsMatch(xqm1, mol3));
            Assert.True(FingerprintsMatch(xqm3, mol3));
        }

        [Fact]
        public void TestBundle1()
        {
            var block = @"
  ChemDraw01052411352D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.247500 0.428587 0.000000 0
M  V30 2 C -0.466950 0.016087 0.000000 0
M  V30 3 C -0.466950 -0.808913 0.000000 0
M  V30 4 C 0.247500 -1.221413 0.000000 0
M  V30 5 N 0.961950 -0.808913 0.000000 0
M  V30 6 C 0.961950 0.016087 0.000000 0
M  V30 7 C -1.181400 0.428587 0.000000 0
M  V30 8 * 0.604725 0.222337 0.000000 0
M  V30 9 Cl 1.181400 1.221413 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 2 7
M  V30 8 1 8 9 ENDPTS=(2 6 1) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
";
            var query = RWMol.MolFromMolBlock(block);
            var xqm = RDKFuncs.createExtendedQueryMol(query, true, true, false);
            var mol = RWMol.MolFromSmiles("Cc1c(cc([nH]1)C(=O)NC2CCN(CC2)c3cc(cc(n3)Cl)C(=O)N)Br");
            var matches = xqm.getSubstructMatches(mol);
            Assert.True(matches.Count > 0);
        }

        [Fact]
        public void TestKekulizeError()
        {
            var block = @"
  ChemDraw01052410132D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.247500 0.428587 0.000000 0
M  V30 2 C -0.466950 0.016087 0.000000 0
M  V30 3 C -0.466950 -0.808913 0.000000 0
M  V30 4 C 0.247500 -1.221413 0.000000 0
M  V30 5 N 0.961950 -0.808913 0.000000 0
M  V30 6 C 0.961950 0.016087 0.000000 0
M  V30 7 C -1.181400 0.428587 0.000000 0
M  V30 8 * 0.604725 0.222337 0.000000 0
M  V30 9 Cl 1.181400 1.221413 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 2 7
M  V30 8 1 8 9 ENDPTS=(2 6 1) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
";
            var query = RWMol.MolFromMolBlock(block);
            ExtendedQueryMol xqm = null;
            var adjustQueryParameters = new AdjustQueryParameters()
            {
                aromatizeIfPossible = false,
                adjustDegree = false,
                adjustSingleBondsBetweenAromaticAtoms = true,
                adjustSingleBondsToDegreeOneNeighbors = true,
                makeDummiesQueries = true,
                adjustConjugatedFiveRings = true
            };
            try
            {
                xqm = RDKFuncs.createExtendedQueryMol(query, true, true, true,
                    adjustQueryParameters);
            }
            catch (ApplicationException ex) when (ex.Message.Contains("MolSanitizeException: Can't kekulize mol"))
            {
                Assert.True(true);
            }

            Assert.Null(xqm);
        }
    }
}