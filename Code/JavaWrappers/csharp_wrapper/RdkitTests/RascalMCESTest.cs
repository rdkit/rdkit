//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using GraphMolWrap;
using Xunit;

namespace RdkitTests;

public class RascalMCESTest
{
    [Fact]
    public void TestTestosteroneVsEstradiol()
    {
        var m1 = RWMol.MolFromSmiles("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C");
        Assert.NotNull(m1);
        var m2 = RWMol.MolFromSmiles("CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O");
        Assert.NotNull(m2);

        var options = new RascalOptions();
        options.similarityThreshold = 0.6;
        var results = RDKFuncs.rascalMCES(m1, m2, options);
        Assert.Equal(1, results.Count);
        var result = results.First();
        var expectedBondMatches = new List<Tuple<int, int>>
        {
            (0, 0).ToTuple(),
            (1, 1).ToTuple(),
            (2, 2).ToTuple(),
            (3, 3).ToTuple(),
            (4, 4).ToTuple(),
            (5, 5).ToTuple(),
            (6, 6).ToTuple(),
            (7, 7).ToTuple(),
            (8, 8).ToTuple(),
            (9, 9).ToTuple(),
            (10, 10).ToTuple(),
            (11, 11).ToTuple(),
            (12, 12).ToTuple(),
            (20, 19).ToTuple(),
            (21, 20).ToTuple(),
            (22, 21).ToTuple()
        };
        var bondMatches = result.getBondMatches().ToList();
        Assert.Equal(16, bondMatches.Count);
        for (int i = 0; i < 16; i++)
        {
            Assert.Equal(expectedBondMatches[i].Item1, bondMatches[i].first);
            Assert.Equal(expectedBondMatches[i].Item2, bondMatches[i].second);
        }

        Assert.Equal(0.4966, result.getSimilarity(), 4);
        var queryMol = RWMol.MolFromSmarts(result.getSmarts());
        Assert.True(m1.hasSubstructMatch(queryMol));
        Assert.True(m2.hasSubstructMatch(queryMol));
    }

    [Fact(Skip = "Works but takes a long time")]
    public void TestRascalCluster()
    {
        var fileName =
            Path.Combine(Environment.GetEnvironmentVariable("RDBASE"),
                "Code", "GraphMol", "RascalMCES", "data", "chembl_1907596.smi");
        var supplier = new SmilesMolSupplier(fileName, "\t", 1, 0, false);
        var molecules = new ROMol_Vect();
        while (!supplier.atEnd())
        {
            molecules.Add(supplier.next());
        }

        var clusterOptions = new RascalClusterOptions();
        clusterOptions.similarityCutoff = 0.7;
        var clusters = RascalApp.RascalCluster(molecules, clusterOptions);

        Assert.Equal(21, clusters.Count);
        var expectedClusterSizes = new[]
        {
            342, 71, 64, 33, 23, 11, 10, 6, 6, 5, 5,
            4, 3, 3, 3, 3, 3, 2, 2, 2, 14
        };
        for (int i = 0; i < 21; i++)
        {
            Assert.Equal(expectedClusterSizes[i], clusters[i].Count);
        }
    }

    [Fact]
    public void TestSmallButina()
    {
        var fileName =
            Path.Combine(Environment.GetEnvironmentVariable("RDBASE"), "Contrib", "Fastcluster", "cdk2.smi");
        var supplier = new SmilesMolSupplier(fileName, "\t", 1, 0, false);
        var molecules = new ROMol_Vect();
        while (!supplier.atEnd())
        {
            molecules.Add(supplier.next());
        }

        var clusters = RascalApp.RascalButinaCluster(molecules);

        Assert.Equal(29, clusters.Count);
        var expectedClusterSizes = new[]
        {
            6, 6, 6, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
        };
        for (int i = 0; i < 29; i++)
        {
            Assert.Equal(expectedClusterSizes[i], clusters[i].Count);
        }
    }

    [Fact]
    public void TestSmall()
    {
        var fileName =
            Path.Combine(Environment.GetEnvironmentVariable("RDBASE"), "Contrib", "Fastcluster", "cdk2.smi");
        var supplier = new SmilesMolSupplier(fileName, "\t", 1, 0, false);
        var molecules = new ROMol_Vect();
        while (!supplier.atEnd())
        {
            molecules.Add(supplier.next());
        }

        var clusters = RascalApp.RascalCluster(molecules);
        Assert.Equal(8, clusters.Count);
        var expectedClusterSizes = new[] { 7, 7, 6, 2, 2, 2, 2, 20 };
        for (int i = 0; i < 8; i++)
        {
            Assert.Equal(expectedClusterSizes[i], clusters[i].Count);
        }
    }
}