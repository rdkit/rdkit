using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GraphMolWrap;
using Xunit;

namespace RdkitTests
{
    public class TestTautomer
    {
        [Fact]
        public void TestTautomerEnumeration()
        {
            var smiles1 = "CCCNC(=N)N";
            var mol1 = RDKFuncs.SmilesToMol(smiles1);
            var smiles2 = "CCCN=C(N)N";
            var mol2 = RDKFuncs.SmilesToMol(smiles2);
            Assert.Equal(0, mol1.getSubstructMatch(mol2).Count);

            var enumerator = new TautomerEnumerator();
            var tautomers1 = enumerator.enumerate(mol1);
            var tautomers2 = enumerator.enumerate(mol2);
            Assert.Equal(2, Convert.ToInt32(tautomers1.size()));
            Assert.Equal(2, Convert.ToInt32(tautomers2.size()));

            var mol1Matches = 0;
            var mol2Matches = 0;
            for (uint i = 0; i < 2; i++)
            {
                var mol1Tautomer = tautomers1.at(i);
                if (mol2.getSubstructMatch(mol1Tautomer).Count > 0) mol2Matches++;
                var mol2Tautomer = tautomers2.at(i);
                if (mol1.getSubstructMatch(mol2Tautomer).Count > 0) mol1Matches++;
            }

            Assert.Equal(1, mol1Matches);
            Assert.Equal(1, mol2Matches);
        }

        [Fact]
        public void TestTautomerCanonicalization()
        {
            var smi =
                "O=C(N)CCC1NC(=O)C2N(C(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC1=O)CCCNC(=N)N)CO)C(O)C)CC(C)C)CCSC)CCC2";
            var mol = RWMol.MolFromSmiles(smi);
            // Record for R
            var molFile = @"
  ChemDraw06032117432D

  0  0  0     0  0              0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 13 0 0 1
M  V30 BEGIN ATOM
M  V30 1 C 1.111765 -1.636387 0.000000 0
M  V30 2 C 0.389632 -1.237272 0.000000 0
M  V30 3 C 0.374442 -0.412275 0.000000 0
M  V30 4 C -0.347639 -0.013245 0.000000 0
M  V30 5 C -0.362863 0.811612 0.000000 0
M  V30 6 N -0.273074 -1.636387 0.000000 0
M  V30 7 O 1.111765 -2.461396 0.000000 0
M  V30 8 N -1.084995 1.210729 0.000000 0
M  V30 9 C -1.100270 2.035672 0.000000 0
M  V30 10 N -0.393589 2.461396 0.000000 0
M  V30 11 N -1.822349 2.434701 0.000000 0
M  V30 12 R1 -0.995042 -1.237157 0.000000 0
M  V30 13 R2 1.822349 -1.217228 0.000000 0
M  V30 14 R3 0.320882 2.048896 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 6 12
M  V30 2 1 1 13
M  V30 3 2 7 1
M  V30 4 1 1 2
M  V30 5 1 2 6
M  V30 6 1 2 3 CFG=1
M  V30 7 1 3 4
M  V30 8 1 4 5
M  V30 9 1 5 8
M  V30 10 2 8 9
M  V30 11 1 9 11
M  V30 12 1 9 10
M  V30 13 1 10 14
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 2)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
";

            var query = RWMol.MolFromMolBlock(molFile);
            foreach (var atom in query.getAtoms())
            {
                if (atom.getAtomicNum() != 0)
                {
                    continue;
                }

                atom.setIsotope(0U);
                atom.setAtomMapNum(0);
            }

            var queryParameters = AdjustQueryParameters.noAdjustments();
            queryParameters.makeDummiesQueries = true;
            var matchParameters = new SubstructMatchParameters
            {
                useChirality = true, specifiedStereoQueryMatchesUnspecified = true, useEnhancedStereo = true
            };

            var cleanupParameters = new CleanupParameters();
            cleanupParameters.tautomerRemoveBondStereo = false;
            cleanupParameters.tautomerRemoveIsotopicHs = false;
            cleanupParameters.tautomerReassignStereo = false;
            cleanupParameters.tautomerRemoveSp3Stereo = false;
            var canonMol = RDKFuncs.canonicalTautomer(mol, cleanupParameters);
            RDKFuncs.addHs(canonMol);
            var canonQuery = RDKFuncs.canonicalTautomer(query, cleanupParameters);
            RDKFuncs.addHs(canonQuery);
            RDKFuncs.adjustQueryProperties(canonQuery, queryParameters);
            var canonMatches = canonMol.getSubstructMatches(canonQuery, matchParameters);
            var numberCanonHits = canonMatches.Count;
            Assert.Equal(1, numberCanonHits);

            RDKFuncs.addHs(mol);
            RDKFuncs.addHs(query);
            RDKFuncs.adjustQueryProperties(query, queryParameters);
            var matches = mol.getSubstructMatches(query, matchParameters);
            var numberHits = matches.Count;
            Assert.Equal(0, numberHits);
        }
    }
}