using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using GraphMolWrap;
using System.Diagnostics;

namespace RDKitCSharpTest
{
    class Program
    {
        static void Main(string[] args)
        {
            // ----- Object creation -----

            Console.WriteLine("Creating some objects:");

            ROMol m1 = RWMol.MolFromSmiles("c1ccccc1");
            Console.WriteLine(" mol: " + m1 + " " + m1.getNumAtoms());
            ROMol m2 = RWMol.MolFromSmiles("c1ccccn1");
            Console.WriteLine(" smi: " + m1 + " " + m1.MolToSmiles());
            Console.WriteLine(" smi2: " + m2 + " " + m2.MolToSmiles());


            ExplicitBitVect fp1 = RDKFuncs.LayeredFingerprintMol(m1);
            ExplicitBitVect fp2 = RDKFuncs.LayeredFingerprintMol(m2);

            Console.WriteLine(" sim: " + RDKFuncs.TanimotoSimilarityEBV(fp1, fp2));

            //rxnTest();
            //smiTest();
            //morganTest();

            ROMol m3 = RWMol.MolFromSmiles("c1ccccc1");
            uint nAtoms = m3.getNumAtoms(true);

            Console.WriteLine("Bulk memory leak test");
            for (uint i = 0; i < 10000; ++i)
            {
                ROMol m4 = RWMol.MolFromSmiles("Clc1cccc(N2CCN(CCC3CCC(CC3)NC(=O)c3cccs3)CC2)c1Cl");
                if ((i % 1000)==0) Console.WriteLine(" Done: " + i);
                m4.Dispose();
                //GC.Collect();
            }

            Console.WriteLine("Goodbye");
        }
    }
}
