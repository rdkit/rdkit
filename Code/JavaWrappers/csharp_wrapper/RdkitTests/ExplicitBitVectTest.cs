using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GraphMolWrap;
using Xunit;

namespace RdkitTests
{
    public class ExplicitBitVectTest
    {
        [Fact]
        public void TestBitOps()
        {
            var mol1 = RWMol.MolFromSmiles("c1ccccc1");
            var mol2 = RWMol.MolFromSmiles("c1ccccn1");
            var fp1 = RDKFuncs.getMorganFingerprintAsBitVect(mol1, 2, 1024);
            var fp2 = RDKFuncs.getMorganFingerprintAsBitVect(mol2, 2, 1024);

            var fp1OnBits = fp1.getOnBits();
            var fp2OnBits = fp2.getOnBits();

            var nCommon = 0;
            var n1 = fp1OnBits.Count;
            var n2 = fp2OnBits.Count;
            foreach (var bit1 in fp1OnBits)
            {
                foreach (var bit2 in fp2OnBits)
                {
                    if (bit1 == bit2)
                    {
                        nCommon++;
                    }
                }
            }

            ExplicitBitVect andFp = new();
            andFp.copy(fp1);
            Assert.Equal(n1, Convert.ToInt32(andFp.getNumOnBits()));
            andFp.andOperator(fp2);
            Assert.Equal(nCommon, Convert.ToInt32(andFp.getNumOnBits()));

            ExplicitBitVect orFp = new();
            orFp.copy(fp1);
            Assert.Equal(n1, Convert.ToInt32(orFp.getNumOnBits()));
            orFp.orOperator(fp2);
            Assert.Equal(n1 + n2 - nCommon, Convert.ToInt32(orFp.getNumOnBits()));

            ExplicitBitVect xorFp = new();
            xorFp.copy(fp1);
            Assert.Equal(n1, Convert.ToInt32(xorFp.getNumOnBits()));
            xorFp.xorOperator(fp2);
            Assert.Equal(n1 + n2 - nCommon*2, Convert.ToInt32(xorFp.getNumOnBits()));


            
        }
    }
}
