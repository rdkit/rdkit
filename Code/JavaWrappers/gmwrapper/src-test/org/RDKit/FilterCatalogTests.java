/* 
 *
 *  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
 *  All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met: 
 *
 *     * Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following 
 *       disclaimer in the documentation and/or other materials provided 
 *       with the distribution.
 *     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
 *       nor the names of its contributors may be used to endorse or promote 
 *       products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package org.RDKit;

import static org.junit.Assert.*;
import org.junit.*;
    
public class FilterCatalogTests extends GraphMolTest {
	private ROMol mol;

	@Before public void setUp() {
	}

	@Test 
	public void test1Basics() {
            FilterCatalog catalog = new FilterCatalog();
            assertEquals(0, catalog.getNumEntries());
            
            //assertEquals(16, catalog.getNumEntries());

	}

        @Test
	public void test2Basics() {

            FilterCatalog catalog = new FilterCatalog(
                     FilterCatalogParams.FilterCatalogs.PAINS_A);
            assertEquals(16, catalog.getNumEntries());

            mol = RWMol.MolFromSmiles("O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2");
            FilterCatalogEntry entry = catalog.getFirstMatch(mol);
            Str_Vect props = entry.getPropList();
            
            String ref  = entry.getProp("Reference");
            String source = entry.getProp("Scope");
            assertEquals(ref,
                         "Baell JB, Holloway GA. New Substructure Filters for " +
                         "Removal of Pan Assay Interference Compounds (PAINS) " +
                         "from Screening Libraries and for Their Exclusion in " +
                         "Bioassays. J Med Chem 53 (2010) 2719D40. " +
                         "doi:10.1021/jm901137j.");

            assertEquals(source, "PAINS filters (family A)");
            assertEquals(entry.getDescription(),"hzone_phenol_A(479)");

            // check the getMatches api point
            FilterCatalogEntry_Vect matches = catalog.getMatches(mol);
            assertEquals(1, matches.size());
            
            for (int entryIdx = 0; entryIdx < matches.size(); ++entryIdx) {
                entry = matches.get(entryIdx);
                String refa  = entry.getProp("Reference");
                String sourcea = entry.getProp("Scope");
                assertEquals(refa,
                             "Baell JB, Holloway GA. New Substructure Filters for " +
                             "Removal of Pan Assay Interference Compounds (PAINS) " +
                             "from Screening Libraries and for Their Exclusion in " +
                             "Bioassays. J Med Chem 53 (2010) 2719D40. " +
                             "doi:10.1021/jm901137j.");
                assertEquals(source, "PAINS filters (family A)");
                assertEquals(entry.getDescription(),"hzone_phenol_A(479)");

                FilterMatch_Vect fmatches = entry.getFilterMatches(mol);
                assertEquals(1, fmatches.size());
                Match_Vect mv = fmatches.get(0).getAtomMatches();
                
                assertEquals(0, mv.get(0).getFirst());
                assertEquals(23, mv.get(0).getSecond());
                assertEquals(1, mv.get(1).getFirst());
                assertEquals(22, mv.get(1).getSecond());
                assertEquals(2, mv.get(2).getFirst());
                assertEquals(20, mv.get(2).getSecond());
                assertEquals(3, mv.get(3).getFirst());
                assertEquals(19, mv.get(3).getSecond());
                assertEquals(4, mv.get(4).getFirst());
                assertEquals(25, mv.get(4).getSecond());
                assertEquals(5, mv.get(5).getFirst());
                assertEquals(24, mv.get(5).getSecond());
                assertEquals(6, mv.get(6).getFirst());
                assertEquals(18, mv.get(6).getSecond());
                assertEquals(7, mv.get(7).getFirst());
                assertEquals(17, mv.get(7).getSecond());
                assertEquals(8, mv.get(8).getFirst());
                assertEquals(16, mv.get(8).getSecond());
                assertEquals(9, mv.get(9).getFirst());
                assertEquals(21, mv.get(9).getSecond());
            }
            
            if (catalog.canSerialize()) {
                byte pickle[] = catalog.Serialize();
                System.out.println(pickle);
                assertTrue(pickle != null);
                FilterCatalog catalog2 = FilterCatalog.Deserialize(pickle);
                assertFalse(catalog2 == null);
                assertEquals(16, catalog2.getNumEntries());
                entry = catalog2.getFirstMatch(mol);
                assertEquals(entry.getDescription(),"hzone_phenol_A(479)");
            }
	}

        @Test
    	public void testRemoveEntry() {
            FilterCatalog catalog = new FilterCatalog(
                FilterCatalogParams.FilterCatalogs.PAINS_A);
            assertEquals(16, catalog.getNumEntries());
            
            mol = RWMol.MolFromSmiles("O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2");
            FilterCatalogEntry entry = catalog.getFirstMatch(mol);
            assertEquals(entry.getDescription(),"hzone_phenol_A(479)");

            assertTrue(catalog.removeEntry(entry));
            FilterCatalogEntry entry_removed = catalog.getFirstMatch(mol);
            assertEquals(entry_removed, null);

            catalog.addEntry(entry);
            entry = catalog.getFirstMatch(mol);
            assertEquals(entry.getDescription(),"hzone_phenol_A(479)");
            
        }

        @Test
	public void testFilterCatalogRunner() {
	    FilterCatalog catalog = new FilterCatalog(
                FilterCatalogParams.FilterCatalogs.PAINS_A);
            assertEquals(16, catalog.getNumEntries());
	    Str_Vect smiles = new Str_Vect(1);
	    smiles.set(0, "O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2");

	    FilterCatalogEntry_VectVect result = RDKFuncs.RunFilterCatalog(catalog, smiles);

	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.FilterCatalogTests");
	}

}
