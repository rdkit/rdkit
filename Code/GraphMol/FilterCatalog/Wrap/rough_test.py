# $Id$
#
#  Copyright (C) 2015  Novartis Institute of BioMedical Research
#         All Rights Reserved
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

""" This is a rough coverage test of the python wrapper for FilterCatalogs

it is intended to be shallow but broad.
"""

from __future__ import print_function
import unittest,os
from rdkit.six.moves import cPickle
from rdkit import RDConfig
from rdkit.RDLogger import logger
logger=logger()
from rdkit import Chem
from rdkit.Chem import FilterCatalog, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalogParams
from rdkit.Chem.FilterCatalog import FilterMatchOps
from rdkit import DataStructs

class TestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test0FilterCatalogEntry(self):
        matcher = FilterCatalog.SmartsMatcher("Aromatic carbon chain")
        self.assertTrue(not matcher.IsValid())

        pat = Chem.MolFromSmarts("c:c:c:c:c")
        matcher.SetPattern(pat)
        matcher.SetMinCount(1)
        
        entry = FilterCatalog.FilterCatalogEntry("Bar", matcher)
        if FilterCatalog.FilterCatalogCanSerialize():
            pickle = entry.Serialize()
        else:
            pickle = None


        self.assertTrue(entry.GetDescription() == "Bar")
        self.assertTrue(matcher.GetMinCount() == 1)
        self.assertTrue(matcher.GetMaxCount() == 2**32-1)
        self.assertTrue(matcher.IsValid())

        entry.SetDescription("Foo")
        self.assertTrue(entry.GetDescription() == "Foo")

        mol = Chem.MolFromSmiles("c1ccccc1")
        self.assertTrue(matcher.HasMatch(mol))


        matcher = FilterCatalog.SmartsMatcher(pat)
        self.assertEqual(str(matcher), "Unnamed SmartsMatcher")
        self.assertTrue(matcher.GetMinCount() == 1)
        self.assertTrue(matcher.HasMatch(mol))
        matches = matcher.GetMatches(mol)

        matcher = FilterCatalog.ExclusionList()
        matcher.SetExclusionPatterns([matcher])
        self.assertTrue(not matcher.HasMatch(mol))
        
        #pat = Chem.MolFromSmarts("c:c:c:c:c")
        #entry.SetOnPattern(pat)
        #entry.SetOffPatterns([pat,pat,pat])
        #self.assertTrue(not entry.HasMatch(pat))

    def test1FilterMatchOps(self):
        mol = Chem.MolFromSmiles("c1ccccc1")

        pat = Chem.MolFromSmarts("c:c:c:c:c")
        matcher = FilterCatalog.SmartsMatcher("Five aromatic carbons", pat)
        self.assertTrue(matcher.GetMinCount() == 1)
        self.assertTrue(matcher.HasMatch(mol))
        matches = matcher.GetMatches(mol)

        matcher2 = FilterCatalog.ExclusionList()
        matcher2.SetExclusionPatterns([matcher])
        self.assertTrue(not matcher2.HasMatch(mol))


        and_match = FilterMatchOps.And(matcher, matcher2)
        self.assertTrue(not and_match.HasMatch(mol))
        not_match = FilterMatchOps.Not(and_match)
        self.assertTrue(not_match.HasMatch(mol))
        or_match = FilterMatchOps.Or(matcher, matcher2)
        self.assertTrue(or_match.HasMatch(mol))

        print(and_match)
        print(or_match)
        print(not_match)
        

    def test2FilterCatalogTest(self):
        tests = ((FilterCatalogParams.FilterCatalogs.PAINS_A, 16),
                 (FilterCatalogParams.FilterCatalogs.PAINS_B, 55),
                 (FilterCatalogParams.FilterCatalogs.PAINS_C, 409),
                 (FilterCatalogParams.FilterCatalogs.PAINS, 409+16+55))

        
        for catalog_idx, num in tests:
            params = FilterCatalog.FilterCatalogParams()
            print ("*"*44)
            print ("Testing:", catalog_idx, int(catalog_idx))
            self.assertTrue(params.AddCatalog(catalog_idx))
            catalog1 = FilterCatalog.FilterCatalog(params)
            
            if FilterCatalog.FilterCatalogCanSerialize():
                pickle = catalog1.Serialize();
                catalog2 = FilterCatalog.FilterCatalog(pickle)
                catalogs = [catalog1, catalog2]
            else:
                catalogs = [catalog1]

            catalogs.append( FilterCatalog.FilterCatalog(catalog_idx) )
            for index,catalog in enumerate(catalogs):
                self.assertEqual(catalog.GetNumEntries(),num)

                if catalog_idx in [FilterCatalogParams.FilterCatalogs.PAINS_A,
                                   FilterCatalogParams.FilterCatalogs.PAINS]:
                    # http://chemistrycompass.com/chemsearch/58909/
                    mol = Chem.MolFromSmiles("O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2")
                    entry = catalog.GetFirstMatch(mol)
                    for key in entry.GetPropList():
                        if key == "Reference":
                            self.assertEquals(entry.GetProp(key),
                                              "Baell JB, Holloway GA. New Substructure Filters for "
                                              "Removal of Pan Assay Interference Compounds (PAINS) "
                                              "from Screening Libraries and for Their Exclusion in "
                                              "Bioassays. J Med Chem 53 (2010) 2719D40. "
                                              "doi:10.1021/jm901137j.")
                        elif key == "Scope":
                            self.assertEquals(entry.GetProp(key), "PAINS filters (family A)")
                            
                    self.assertEqual(entry.GetDescription(),"hzone_phenol_A(479)")
                    result = catalog.GetMatches(mol)
                    self.assertEquals(len(result), 1)
                    
                    for entry in result:
                        for filtermatch in entry.GetFilterMatches(mol):
                            self.assertEquals(str(filtermatch.filterMatch), "hzone_phenol_A(479)")
                            atomPairs = [tuple(x) for x in filtermatch.atomPairs]
                            self.assertEquals(atomPairs,
                                              [(0, 23), (1, 22),
                                               (2, 20), (3, 19),
                                               (4, 25), (5, 24),
                                               (6, 18), (7, 17),
                                               (8, 16), (9, 21)])
                            
                elif catalog_idx == FilterCatalogParams.FilterCatalogs.PAINS_B:
                    mol = Chem.MolFromSmiles("FC(F)(F)Oc1ccc(NN=C(C#N)C#N)cc1") # CHEMBL457504
                    entry = catalog.GetFirstMatch(mol)
                    self.assertTrue(entry)
                    self.assertEquals(entry.GetDescription(), "cyano_imine_B(17)")

                elif catalog_idx == FilterCatalogParams.FilterCatalogs.PAINS_C:
                    mol = Chem.MolFromSmiles("O=C1C2OC2C(=O)c3cc4CCCCc4cc13") # CHEMBL476649
                    entry = catalog.GetFirstMatch(mol)
                    self.assertTrue(entry)
                    self.assertEquals(entry.GetDescription(), "keto_keto_gamma(5)")


    def test3ExclusionFilter(self):
        mol = Chem.MolFromSmiles("c1ccccc1")

        pat = Chem.MolFromSmarts("c:c:c:c:c")
        matcher = FilterCatalog.SmartsMatcher("Five aromatic carbons", pat)
        self.assertTrue(matcher.GetMinCount() == 1)
        self.assertTrue(matcher.HasMatch(mol))
        matches = matcher.GetMatches(mol)

        exclusionFilter = FilterCatalog.ExclusionList()
        exclusionFilter.AddPattern(matcher)
        self.assertFalse(exclusionFilter.HasMatch(mol))
        
        matches2 = exclusionFilter.GetMatches(mol)

        self.assertTrue(matches)
        self.assertFalse(matches2)

    def test4CountTests(self):
        matcher = FilterCatalog.SmartsMatcher("Carbon", "[#6]", 0, 2)
        m = Chem.MolFromSmiles("N")
        self.assertTrue(matcher.HasMatch(m))
        m = Chem.MolFromSmiles("C")
        self.assertTrue(matcher.HasMatch(m))
        m = Chem.MolFromSmiles("CC")
        self.assertTrue(matcher.HasMatch(m))
        m = Chem.MolFromSmiles("CCC")
        self.assertFalse(matcher.HasMatch(m))

        matcher = FilterCatalog.SmartsMatcher("Carbon", "[#6]", 1, 2)
        m = Chem.MolFromSmiles("N")
        self.assertFalse(matcher.HasMatch(m))

        
    def testZinc(self):
        params = FilterCatalog.FilterCatalogParams(FilterCatalogParams.FilterCatalogs.ZINC)
        catalog = FilterCatalog.FilterCatalog(params)
        self.assertTrue(catalog.GetNumEntries())

        m = Chem.MolFromSmiles("C"*41)
        entry = catalog.GetFirstMatch(m)
        self.assertTrue(entry.GetDescription(), "Non-Hydrogen_atoms")

        m = Chem.MolFromSmiles("CN"*20)
        entry = catalog.GetFirstMatch(m)
        self.assertEquals(catalog.GetFirstMatch(m), None)
        
    def testSmartsMatcherAPI(self):
        sm = FilterCatalog.SmartsMatcher("Too many carbons", "[#6]", 40+1)
        sm2 = FilterCatalog.SmartsMatcher("ok # carbons", "[#6]", 0, 40)
        sm3 = FilterCatalog.FilterMatchOps.Not(sm2)
        
        m = Chem.MolFromSmiles("C"*40)
        self.assertFalse(sm.HasMatch(m))
        self.assertTrue(sm2.HasMatch(m))
        self.assertFalse(sm3.HasMatch(m))
        
        m = Chem.MolFromSmiles("C"*41)
        self.assertTrue(sm.HasMatch(m))
        self.assertFalse(sm2.HasMatch(m))
        self.assertTrue(sm3.HasMatch(m))

    def testAddEntry(self):
        sm = FilterCatalog.SmartsMatcher("Too many carbons", "[#6]", 40+1)
        entry = FilterCatalog.FilterCatalogEntry("Bar", sm)
        fc = FilterCatalog.FilterCatalog()
        fc.AddEntry(entry)
        del entry
        del fc

    def testRemoveEntry(self):
        params = FilterCatalog.FilterCatalogParams(FilterCatalogParams.FilterCatalogs.ZINC)
        catalog = FilterCatalog.FilterCatalog(params)
        entry = catalog.GetEntryWithIdx(10)
        desc = entry.GetDescription()
        count = 0
        descs = set( [ catalog.GetEntryWithIdx(i).GetDescription()
                       for i in range (catalog.GetNumEntries())])
        for i in range(catalog.GetNumEntries()):
            if catalog.GetEntryWithIdx(i).GetDescription() == desc:
                count += 1
        print ("Count", count)
        sz = catalog.GetNumEntries()
        print ("*"*44)
        print (dir(catalog))
        self.assertTrue(catalog.RemoveEntry(entry))
        del entry
        self.assertTrue(catalog.GetNumEntries() == sz-1)

        descs2 = set( [ catalog.GetEntryWithIdx(i).GetDescription()
                        for i in range(catalog.GetNumEntries())])
        print (descs - descs2)

        newcount = 0
        for i in range(catalog.GetNumEntries()):
            if catalog.GetEntryWithIdx(i).GetDescription() == desc:
                newcount += 1
        self.assertEquals(count, newcount+1)
        

    def testPyFilter(self):
        class MyFilterMatcher(FilterCatalog.FilterMatcher):
            def IsValid(self):
                return True

            def HasMatch(self, mol):
                return True
            
            def GetMatches(self, mol, vect):
                v = FilterCatalog.MatchTypeVect()
                v.append(FilterCatalog.IntPair(1,1))
                match = FilterCatalog.FilterMatch(self,v)
                vect.append(match)
                return True

        func = MyFilterMatcher("FilterMatcher")
        self.assertEquals(func.GetName(), "FilterMatcher")
        mol = Chem.MolFromSmiles("c1ccccc1")
        self.assertEquals(func.HasMatch(mol), True)

        or_match = FilterMatchOps.Or(func, func)
        self.assertEquals( [[tuple(x) for x in filtermatch.atomPairs]
                          for filtermatch in or_match.GetMatches(mol)],
                           [[(1,1)],[(1,1)]])


        not_match = FilterMatchOps.Not(func)
        print(not_match)
        self.assertEquals(not_match.HasMatch(mol), False)
        # test memory
        del func
        
        self.assertEquals(not_match.HasMatch(mol), False)
        self.assertEquals( [[tuple(x) for x in filtermatch.atomPairs]
                          for filtermatch in not_match.GetMatches(mol)],
                           [])


        entry = FilterCatalog.FilterCatalogEntry("Bar", MyFilterMatcher("FilterMatcher"))
        fc = FilterCatalog.FilterCatalog()
        fc.AddEntry(entry)

        catalogEntry = fc.GetFirstMatch(mol)
        print (catalogEntry.GetDescription())

    def testMWFilter(self):
        class MWFilter(FilterCatalog.FilterMatcher):
            def __init__(self, minMw, maxMw):
                FilterCatalog.FilterMatcher.__init__(self, "MW violation")
                self.minMw = minMw
                self.maxMw = maxMw

            def IsValid(self):
                return True

            def HasMatch(self, mol):
                mw = rdMolDescriptors.CalcExactMolWt(mol)
                return not self.minMw <= mw <= self.maxMw

        entry = FilterCatalog.FilterCatalogEntry("MW Violation", MWFilter(100,500))
        fc = FilterCatalog.FilterCatalog()
        fc.AddEntry(entry)
        self.assertTrue(entry.GetDescription() == "MW Violation")

        mol = Chem.MolFromSmiles("c1ccccc1")
        catalogEntry = fc.GetFirstMatch(mol)



        
if __name__ == '__main__':
    unittest.main()
