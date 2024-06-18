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
    }
}