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
    
	public void test2Basics() {

            FilterCatalog catalog = new FilterCatalog(
                     FilterCatalogParams.FilterCatalogs.PAINS_A);
            assertEquals(16, catalog.getNumEntries());

            mol = RWMol.MolFromSmiles("O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2");
            FilterCatalogEntry entry = catalog.getFirstMatch(mol);
            //Str_Vect props = entry.getPropList();
            
            //for (int i=0; i< props.size(); ++i) {
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

            /*           
            if (catalog.CanSerialize()) {
                String pickle = catalog.Serialize();
                FilterCatalog catalog2 = FilterCatalog.FilterCatalog(pickle);
                assertEquals(16, catalog2.getNumEntries());
            }
            */

	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.FilterCatalogTests");
	}

}
