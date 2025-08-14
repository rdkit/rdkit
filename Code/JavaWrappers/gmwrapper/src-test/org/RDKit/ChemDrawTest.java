/* 
 * $Id: ChemTests.java 131 2011-01-20 22:01:29Z ebakke $
 *
 *  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
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
import java.io.File;

import org.junit.Test;

public class ChemDrawTest extends GraphMolTest {
    @Test public void testChemDrawReader() {
	String rdpath = System.getenv("RDBASE");
	if (rdpath == null)
	    org.junit.Assert.fail("No definition for RDBASE");
	File base = new File(rdpath);
	File testFile = new File(base, "Code" + File.separator + "GraphMol"
				 + File.separator + "test_data" + File.separator +
				 "CDXML" + File.separator + "beta-cypermethrin.cdxml");
	String fn = testFile.getAbsolutePath();
	RWMol_Vect prods = RDKFuncs.MolsFromChemDrawFile(fn);
	assertEquals(prods.size(), 1);
	for(int idx = 0; idx < prods.size(); idx++) {
	    if(idx == 0) {
            System.out.print(prods.get(idx).MolToSmiles(true));
            System.out.print("\n");
		assertEquals(prods.get(idx).MolToSmiles(true), "CC1(C)[C@H](C=C(Cl)Cl)[C@H]1C(=O)O[C@@H](C#N)c1cccc(Oc2ccccc2)c1");
	    }
	}
	ChemDrawParserParams params = new ChemDrawParserParams();
	prods = RDKFuncs.MolsFromChemDrawFile(fn, params);
	assertEquals(prods.size(), 1);
	for(int idx = 0; idx < prods.size(); idx++) {
	    if(idx == 0) {
            System.out.print(prods.get(idx).MolToSmiles(true));
            System.out.print("\n");
		assertEquals(prods.get(idx).MolToSmiles(true), "CC1(C)[C@H](C=C(Cl)Cl)[C@H]1C(=O)O[C@@H](C#N)c1cccc(Oc2ccccc2)c1");
	    }
	}

	testFile = new File(base, "Code" + File.separator + "GraphMol"
				 + File.separator + "test_data" + File.separator +
				 "CDXML" + File.separator + "ring-stereo1.cdx");
	fn = testFile.getAbsolutePath();
	params = new ChemDrawParserParams(true, true, CDXFormat.CDX);
	prods = RDKFuncs.MolsFromChemDrawFile(fn, params);
	assertEquals(prods.size(), 1);

	params = new ChemDrawParserParams(true, true, CDXFormat.CDXML);
	boolean e = false;
	try {
	    prods = RDKFuncs.MolsFromChemDrawFile(fn, params);
	} catch(GenericRDKitException ex) {
	    e = true;
	}
	assertEquals(true, e);
	
    }    


    public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("org.RDKit.ChemDrawTest");
    }

}
