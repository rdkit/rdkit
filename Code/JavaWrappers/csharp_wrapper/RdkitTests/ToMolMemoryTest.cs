using System;
using System.Diagnostics;
using GraphMolWrap;
using Xunit;

namespace RdkitTests
{
    public class ToMolMemoryTest
    {
        private static readonly long OneHundredMB = 1024 * 1024 * 100;
        private static readonly long TwoHundredMB = OneHundredMB * 2;

        private static void gc()
        {
            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();
        }


        [Fact]
        public void TestSmilesToMolMemoryUsage()
        {
            string smi =
                "CC(C)C[C@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CO)NC(=O)[C@@H]1CCCN1C(=O)[C@H](CCC(N)=O)NC(=O)[C@@H](N)CS)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CS)C(=O)O";

            var before = Process.GetCurrentProcess().VirtualMemorySize64;
            var privateBefore = Process.GetCurrentProcess().PrivateMemorySize64;
            long after;
            long privateAfter;
            for (int i = 0; i < 500; ++i)
            {
                RWMol mol = RDKFuncs.SmilesToMol(smi);
                RWMol mol2 = RWMol.MolFromSmiles(smi);
                if (i % 50 == 0)
                {
                    gc();
                    after = Process.GetCurrentProcess().VirtualMemorySize64;
                    Assert.True(after - before < TwoHundredMB);
                    privateAfter = Process.GetCurrentProcess().PrivateMemorySize64;
                    Assert.True(privateAfter - privateBefore < OneHundredMB);
                }
            }

            gc();
            after = Process.GetCurrentProcess().VirtualMemorySize64;
            Assert.True(after - before < TwoHundredMB);
            privateAfter = Process.GetCurrentProcess().PrivateMemorySize64;
            Assert.True(privateAfter - privateBefore < OneHundredMB);
        }

        [Fact]
        public void TestMolBlockToMolMemoryUsage()
        {
            var block = @"
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
";
            var before = Process.GetCurrentProcess().VirtualMemorySize64;
            long after;
            for (int i = 0; i < 500; ++i)
            {
                RWMol mol = RDKFuncs.MolBlockToMol(block);
                if (i % 50 == 0)
                {
                    gc();
                    after = Process.GetCurrentProcess().VirtualMemorySize64;
                    Assert.True(after - before < TwoHundredMB);
                }
            }

            gc();
            after = System.Diagnostics.Process.GetCurrentProcess().PeakVirtualMemorySize64;
            Assert.True(after - before < TwoHundredMB);
        }
    }
}