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
    }
}