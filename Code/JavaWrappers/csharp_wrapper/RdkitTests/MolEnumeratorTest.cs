using GraphMolWrap;
using Xunit;

namespace RdkitTests
{

    public class MolEnumeratorTest
    {
        [Fact]
        public void TestEnumerator()
        {
            var block = @"
  Mrv2008 07102308312D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 10 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.2917 4.9984 0 0
M  V30 2 C -5.6253 4.2284 0 0
M  V30 3 C -5.6253 2.6883 0 0
M  V30 4 C -4.2917 1.9183 0 0
M  V30 5 C -2.958 2.6883 0 0
M  V30 6 C -2.958 4.2284 0 0
M  V30 7 * -3.6248 4.6134 0 0
M  V30 8 Cl -3.6248 6.9234 0 0
M  V30 9 * -5.6253 3.4583 0 0
M  V30 10 Cl -7.2521 3.4583 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 2 4 5
M  V30 5 1 5 6
M  V30 6 2 1 6
M  V30 7 1 7 8 ENDPTS=(2 1 6) ATTACH=ANY
M  V30 8 1 9 10 ENDPTS=(2 2 3) ATTACH=ANY
M  V30 END BOND
M  V30 END CTAB
M  END
";
            var mol = RWMol.MolFromMolBlock(block);
            var bundle = RDKFuncs.enumerate(mol);
            Assert.Equal(4U, bundle.size());
        }
    }
}