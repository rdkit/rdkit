using System.Collections.Generic;
using GraphMolWrap;
using Xunit;


namespace RdkitTests
{
    public class TestMolStandardize
    {
        [Fact]
        public void TestNormalize()
        {
            var smiles = "C[S+2]([O-])([O-])C([O-])C(=O)O";
            var mol = RWMol.MolFromSmiles(smiles);
            var normalizedMol = RDKFuncs.normalize(mol);
            Assert.Equal("CS(=O)(=O)C([O-])C(=O)O", normalizedMol.MolToSmiles());
            Assert.Equal(smiles, mol.MolToSmiles());
            RDKFuncs.normalizeInPlace(mol);
            Assert.Equal("CS(=O)(=O)C([O-])C(=O)O", mol.MolToSmiles());
        }


        [Fact]
        public void TestCleanup()
        {
            var smiles = "O=C(O)[C@]([O-])(O)Cl";
            var mol = RWMol.MolFromSmiles(smiles);
            var cleanedMol = RDKFuncs.cleanup(mol);
            Assert.Equal("O=C([O-])C(O)(O)Cl", cleanedMol.MolToSmiles());
            Assert.Equal(smiles, mol.MolToSmiles());
            RDKFuncs.cleanupInPlace(mol);
            Assert.Equal("O=C([O-])C(O)(O)Cl", mol.MolToSmiles());
        }


        [Fact]
        public void TestCanonicalTautomer()
        {
            var smiles = "CP(C)O";
            var mol = RWMol.MolFromSmiles(smiles);
            var tautomer = RDKFuncs.canonicalTautomer(mol);
            Assert.Equal("C[PH](C)=O", tautomer.MolToSmiles());
            Assert.Equal(smiles, mol.MolToSmiles());
            RDKFuncs.canonicalTautomerInPlace(mol);
            Assert.Equal("C[PH](C)=O", mol.MolToSmiles());
        }

        // Test that standardization rules for a guanidinium group work as expected when matching an Arginine query to a set of peptides
        [Fact]
        public void TestRMatchNormalize()
        {
            var cleanupParameters = new CleanupParameters();
            cleanupParameters.doCanonical = true;
            var replacements = cleanupParameters.normalizationData;
            var replacement = new String_Pair("Standardize ARG",
                "[#6:1][#6:2][#6:3][#7H1:4]-[#6X3:5](=[#7:6])-[#7:7]>>[#6:1][#6:2][#6:3][#7H0:4]=[#6X3:5](-[#7:6])-[#7:7]");
            replacements.Add(replacement);
            cleanupParameters.normalizationData = replacements;
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
            var q = RWMol.MolFromMolBlock(molFile);
            q = RDKFuncs.normalize(q, cleanupParameters);
            RDKFuncs.addHs(q);
            foreach (var atom in q.getAtoms())
            {
                if (atom.getAtomicNum() != 0) continue;
                atom.setAtomMapNum(0);
                atom.setIsotope(0);
                atom.clearProp("_MolFileRLabel");
            }

            var queryParameters = AdjustQueryParameters.noAdjustments();
            queryParameters.makeDummiesQueries = true;
            RDKFuncs.adjustQueryProperties(q, queryParameters);
            var matchParameters = new SubstructMatchParameters
            {
                useChirality = true, specifiedStereoQueryMatchesUnspecified = true, useEnhancedStereo = true
            };
            List<string> smiles = new();
            smiles.Add(
                "O=C(N)CCC1NC(=O)C2N(C(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC1=O)CCCNC(=N)N)CO)C(O)C)CC(C)C)CCSC)CCC2");
            smiles.Add(
                "C[C@H](N[H])C(=O)N[C@H](C(=O)N[C@@H]1C(=O)N([C@@H](C)C(=O)N[C@H](C)C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C=O)CC(O)=O)CSSC1)[C@@H](C)O)CCCCN[H])[C@@H](C)C(=O)O)CCCNC(N)=N");
            smiles.Add(
                "CSCC[C@@H]1NC(=O)[C@H](CC(C)C)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CO)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@@H]2CCCN2C1=O");
            smiles.Add(
                "C[C@H](N)C(=O)N[C@@H](CCCN=C(N)N)C(=O)N[C@H]1CSSC[C@@H](C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](C)C(=O)O)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](C)NC1=O");
            foreach (var smi in smiles)
            {
                var mol = RWMol.MolFromSmiles(smi);

                mol = RDKFuncs.normalize(mol, cleanupParameters);
                RDKFuncs.addHs(mol);

                var normalizedHits = mol.getSubstructMatches(q, matchParameters);
                Assert.Equal(1, normalizedHits.Count);
            }
        }
    }

}