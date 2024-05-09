// Linux:
// compile with
// mcs -platform:x64 -r:../RDKit2DotNet.dll -out:MolToFromByteArray.exe  MolToFromByteArray.cs
// and run with
// LD_LIBRARY_PATH=..:$RDBASE/lib:$LD_LIBRARY_PATH MONO_PATH=.. mono MolToFromByteArray.exe

using System.IO;
using System.Diagnostics;
using GraphMolWrap;
using Xunit;

namespace RdkitTests
{
    public class MolToFromByteArrayTest
    {
        [Fact]
        public void TestMolToFromByteArray()
        {
            string smi = "CN(C)c1ccc2c(=O)cc[nH]c2c1";
            string pklFileName = "quinolone.pkl";
            {
                ROMol mol = RWMol.MolFromSmiles(smi);
                byte[] pkl = mol.ToByteArray();
                File.WriteAllBytes(pklFileName, pkl);
                mol.Dispose();
            }
            {
                byte[] pkl = File.ReadAllBytes(pklFileName);
                ROMol mol = ROMol.FromByteArray(pkl);
                Assert.Equal(smi, mol.MolToSmiles());
                mol.Dispose();
                File.Delete(pklFileName);
            }
        }
    }
}