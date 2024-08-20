//
//  Copyright (C) 2024 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

using System;
using System.IO;
using System.Linq;
using GraphMolWrap;
using Xunit;

namespace RdkitTests
{
    public class TestSubstanceGroup
    {
        [Fact]
        public void TestReadSRUSGroup()
        {
            var fName =
                Path.Combine(Environment.GetEnvironmentVariable("RDBASE")!,
                    "Code/GraphMol/FileParsers/sgroup_test_data/repeat_groups_query1.mol");
            var block = File.ReadAllText(fName);
            var mol = RWMol.MolFromMolBlock(block);

            var numberGroups = RDKFuncs.getSubstanceGroupCount(mol);
            Assert.Equal(1U, numberGroups);
            var sruGroup = RDKFuncs.getSubstanceGroupWithIdx(mol, 0);
            var atoms = sruGroup.getAtoms().ToArray();
            Assert.Equal(new int[] { 0 }, atoms);
            var bonds = sruGroup.getBonds().ToArray();
            Assert.Equal(new int[] { 0, 5 }, bonds);
            var index = Convert.ToInt32(sruGroup.getUIntProp("index"));
            Assert.Equal(1, index);
            var connect = sruGroup.getStringProp("CONNECT");
            Assert.Equal("HT", connect);
            var label = sruGroup.getStringProp("LABEL");
            Assert.Equal("1-3", label);

            var props = sruGroup.getPropList().ToArray();

            Assert.True(props.Contains("XBHEAD"));
            var xbBonds = sruGroup.getUIntVectProp("XBHEAD").Select(Convert.ToInt32).ToArray();
            Assert.Equal(new int[] { 5, 0 }, xbBonds);

            Assert.True(props.Contains("XBCORR"));
            var xbCorr = sruGroup.getUIntVectProp("XBCORR").Select(Convert.ToInt32).ToArray();
            Assert.Equal(new int[] { 5, 5, 0, 0 }, xbCorr);
        }
    }
}