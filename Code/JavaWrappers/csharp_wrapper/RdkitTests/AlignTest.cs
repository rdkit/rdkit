//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

using System;
using System.IO;
using GraphMolWrap;
using Xunit;

namespace RdkitTests;

// Based on Code/GraphMol/MolAlign/Wrap/testMolAlign.py
public class AlignTest
{
    private static string AlignDataFile(string name)
    {
        return Path.Combine(Environment.GetEnvironmentVariable("RDBASE")!, "Code", "GraphMol", "MolAlign",
            "test_data", name);
    }

    private static void AssertClose(double expected, double actual, double tol)
    {
        Assert.InRange(actual, expected - tol, expected + tol);
    }

    private static Match_Vect OirAtomMap()
    {
        var atomMap = new Match_Vect();
        foreach (var (prb, refIdx) in new[] { (18, 27), (13, 23), (21, 14), (24, 7), (9, 19), (16, 30) })
        {
            atomMap.Add(new Int_Pair(prb, refIdx));
        }
        return atomMap;
    }

    [Fact]
    public void TestAlignMol()
    {
        var mol1 = RWMol.MolFromMolFile(AlignDataFile("1oir.mol"));
        var mol2 = RWMol.MolFromMolFile(AlignDataFile("1oir_conf.mol"));
        var rmsd = RDKFuncs.alignMol(mol2, mol1);
        AssertClose(0.6578, rmsd, 0.0001);

        // the probe should now have the coordinates of the pre-aligned file
        var mol3 = RWMol.MolFromMolFile(AlignDataFile("1oir_trans.mol"));
        var conf2 = mol2.getConformer();
        var conf3 = mol3.getConformer();
        for (uint i = 0; i < mol2.getNumAtoms(); ++i)
        {
            var p2 = conf2.getAtomPos(i);
            var p3 = conf3.getAtomPos(i);
            AssertClose(p3.x, p2.x, 0.001);
            AssertClose(p3.y, p2.y, 0.001);
            AssertClose(p3.z, p2.z, 0.001);
        }

        var trans = new Transform3D();
        rmsd = RDKFuncs.getAlignmentTransform(mol2, mol1, trans);
        AssertClose(0.6579, rmsd, 0.0001);
    }

    [Fact]
    public void TestAlignMolAtomMapAndWeights()
    {
        var mol1 = RWMol.MolFromMolFile(AlignDataFile("1oir.mol"));
        var mol2 = RWMol.MolFromMolFile(AlignDataFile("1oir_conf.mol"));
        var atomMap = OirAtomMap();
        var rmsd = RDKFuncs.alignMol(mol2, mol1, 0, 0, atomMap);
        AssertClose(0.8525, rmsd, 0.0001);

        mol2 = RWMol.MolFromMolFile(AlignDataFile("1oir_conf.mol"));
        var wts = new DoubleVector(6, 1.0);
        wts.setVal(5, 2.0);
        rmsd = RDKFuncs.alignMol(mol2, mol1, 0, 0, atomMap, wts);
        AssertClose(0.9513, rmsd, 0.0001);
    }

    [Fact]
    public void TestAlignMolConformers()
    {
        var mol = RWMol.MolFromSmiles("C1CC1CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4");
        var cids = DistanceGeom.EmbedMultipleConfs(mol, 10, 30, 100);
        Assert.Equal(10, cids.Count);
        var aids = new UInt_Vect();
        foreach (var aid in new uint[] { 12, 13, 14, 15, 16, 17, 18 })
        {
            aids.Add(aid);
        }
        var rmsVals = new Double_Vect();
        RDKFuncs.alignMolConformers(mol, aids, null, null, false, 50, rmsVals);
        Assert.Equal((int)mol.getNumConformers() - 1, rmsVals.Count);
        foreach (var rms in rmsVals)
        {
            Assert.True(rms >= 0.0);
        }
    }

    [Fact]
    public void TestBestRMS()
    {
        var suppl = new SDMolSupplier(AlignDataFile("probe_mol.sdf"), true, false);
        suppl.moveTo(1);
        var prb = suppl.next();
        var refMol = suppl.next();
        var prbCopy = new ROMol(prb);

        var rmsdInPlace = RDKFuncs.CalcRMS(new ROMol(prb), refMol);
        AssertClose(2.6026, rmsdInPlace, 0.001);
        // alignMol() would give 2.50561 here
        var rmsd = RDKFuncs.getBestRMS(prb, refMol);
        AssertClose(2.43449, rmsd, 0.001);

        var bestTrans = new Transform3D();
        var bestMatch = new Match_Vect();
        var rmsdCopy = RDKFuncs.getBestAlignmentTransform(prbCopy, refMol, bestTrans, bestMatch);
        AssertClose(rmsd, rmsdCopy, 0.001);
        Assert.Equal((int)refMol.getNumAtoms(), bestMatch.Count);
    }

    [Fact]
    public void TestBestAlignmentParams()
    {
        var suppl = new SDMolSupplier(AlignDataFile("probe_mol.sdf"), true, false);
        suppl.moveTo(1);
        var prb = suppl.next();
        var refMol = suppl.next();

        // Hs are ignored by default
        var ps = new BestAlignmentParams();
        Assert.True(ps.ignoreHs);
        AssertClose(1.8100, RDKFuncs.getBestRMS(new ROMol(prb), refMol, ps), 0.001);

        ps.ignoreHs = false;
        AssertClose(2.43449, RDKFuncs.getBestRMS(new ROMol(prb), refMol, ps), 0.001);

        var bestTrans = new Transform3D();
        var bestMatch = new Match_Vect();
        var rmsd = RDKFuncs.getBestAlignmentTransform(new ROMol(prb), refMol, bestTrans, bestMatch, ps);
        AssertClose(2.43449, rmsd, 0.001);
        Assert.Equal((int)refMol.getNumAtoms(), bestMatch.Count);
    }

    [Fact]
    public void TestGetAllConformerBestRMSToRef()
    {
        var prbMol = RWMol.MolFromSmiles("OCCCN1CCN(C)CC1");
        DistanceGeom.EmbedMultipleConfs(prbMol, 5, 30, 42);
        Assert.Equal(5u, prbMol.getNumConformers());
        // the reference only has a copy of the first probe conformer
        var refMol = new ROMol(prbMol, false, 0);
        Assert.Equal(1u, refMol.getNumConformers());

        var ps = new BestAlignmentParams();
        var rmsds = RDKFuncs.getAllConformerBestRMSToRef(prbMol, refMol, ps);
        Assert.Equal(5, rmsds.Count);
        AssertClose(0.0, rmsds[0], 0.0001);
        for (var i = 0; i < rmsds.Count; ++i)
        {
            AssertClose(RDKFuncs.getBestRMS(new ROMol(prbMol), refMol, ps, i, 0), rmsds[i], 0.0001);
        }
    }

    [Fact]
    public void TestBestRMSConjugatedGroups()
    {
        var mol = RWMol.MolFromSmiles(
            "CCC(=O)[O-] |(-1.11,0.08,-0.29;0.08,-0.18,0.58;1.34,0.03,-0.16;1.74,1.22,-0.32;2.06,-1.04,-0.66)|");
        var qry = RWMol.MolFromSmiles(
            "CCC([O-])=O |(-1.11,0.08,-0.29;0.08,-0.18,0.58;1.34,0.03,-0.16;1.74,1.22,-0.32;2.06,-1.04,-0.66)|");
        AssertClose(0.0, RDKFuncs.getBestRMS(qry, mol), 0.001);
        AssertClose(0.747, RDKFuncs.getBestRMS(qry, mol, -1, -1, new Match_Vect_Vect(), 1000000, false), 0.001);
    }

    [Fact]
    public void TestGetAllConformerBestRMS()
    {
        var mol = RWMol.MolFromSmiles("OCCCN1CCN(C)CC1");
        DistanceGeom.EmbedMultipleConfs(mol, 5, 30, 42);
        var nconfs = (int)mol.getNumConformers();
        Assert.Equal(5, nconfs);
        var origVals = RDKFuncs.getAllConformerBestRMS(mol);
        Assert.Equal(nconfs * (nconfs - 1) / 2, origVals.Count);

        var newVals = RDKFuncs.getAllConformerBestRMS(mol, 4);
        Assert.Equal(origVals.Count, newVals.Count);
        for (var i = 0; i < origVals.Count; ++i)
        {
            AssertClose(origVals[i], newVals[i], 0.0001);
        }

        AssertClose(origVals[0], RDKFuncs.getBestRMS(mol, mol, 0, 1), 0.0001);
    }

    private static void CheckO3AOnRefE2(O3A.AtomTypeScheme atomTypes, double expectedScore, double expectedRMSD)
    {
        var suppl = new SDMolSupplier(AlignDataFile("ref_e2.sdf"), true, false);
        suppl.moveTo(48);
        var refMol = suppl.next();
        suppl.reset();
        var cumScore = 0.0;
        var cumMsd = 0.0;
        var nMols = 0;
        while (!suppl.atEnd())
        {
            var prbMol = suppl.next();
            using var o3a = new O3A(prbMol, refMol, atomTypes);
            cumScore += o3a.score();
            var rmsd = o3a.align();
            cumMsd += rmsd * rmsd;
            ++nMols;
        }
        cumMsd /= nMols;
        AssertClose(expectedScore, cumScore, 1.0);
        AssertClose(expectedRMSD, Math.Sqrt(cumMsd), 0.001);
    }

    [Fact]
    public void TestO3AMMFF()
    {
        CheckO3AOnRefE2(O3A.AtomTypeScheme.MMFF94, 6942, 0.345);
    }

    [Fact]
    public void TestO3ACrippen()
    {
        CheckO3AOnRefE2(O3A.AtomTypeScheme.CRIPPEN, 4918, 0.304);
    }

    [Fact]
    public void TestO3AMatchesAndTransform()
    {
        var suppl = new SDMolSupplier(AlignDataFile("ref_e2.sdf"), true, false);
        var refMol = suppl.next();
        var prbMol = suppl.next();
        using var o3a = new O3A(prbMol, refMol);
        var matches = o3a.matches();
        Assert.True(matches.Count > 0);
        Assert.Equal((uint)matches.Count, o3a.weights().size());
        var trans = new Transform3D();
        var transRmsd = o3a.trans(trans);
        AssertClose(transRmsd, o3a.align(), 1e-6);
    }

    [Fact]
    public void TestO3AMissingMMFFParams()
    {
        var m1 = RWMol.MolFromSmiles("c1ccccc1Cl");
        DistanceGeom.EmbedMolecule(m1);
        var m2 = RWMol.MolFromSmiles("c1ccccc1B(O)O");
        DistanceGeom.EmbedMolecule(m2);
        Assert.ThrowsAny<Exception>(() => new O3A(m1, m2));
    }
}
