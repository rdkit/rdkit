/* 
 * $Id: SuppliersTests.java 131 2011-01-20 22:01:29Z ebakke $
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

import java.io.*;
import java.util.ArrayList;

import org.junit.*;

public class SuppliersTests extends GraphMolTest {

	private ArrayList<String> tmpFiles = new ArrayList<String>();
	private File baseTestPath;

	@Before 
	public void setUp() {
		File base = getRdBase();
		baseTestPath = new File(base, "rdkit" + File.separator + "VLib" + File.separator + 
				"NodeLib" + File.separator + "test_data");

	}

	@After
	public void tearDown() {
		for (String fn : tmpFiles) {
			new File(fn).delete();
		}
	}

	@Test
	public void test1SDSupplier() {
		File fileN = new File(baseTestPath, "NCI_aids.10.sdf");
		SDMolSupplier suppl = new SDMolSupplier(fileN.getPath());
		ArrayList<ROMol> ms = new ArrayList<ROMol>();
		ROMol m;
		do { 
			m = suppl.next();
			ms.add(m);
		}
		while (!suppl.atEnd());
		assertEquals( 10,ms.size() );
		ms.clear();

		suppl.reset();
		ms.clear();
		do { 
			m = suppl.next();
			ms.add(m);
		}
		while (!suppl.atEnd());
		assertEquals( 10,ms.size() );

		suppl.reset();
		m = suppl.next();
		assertEquals("48",m.getProp("_Name"));
		assertEquals("48",m.getProp("NSC"));
		assertEquals("15716-70-8",m.getProp("CAS_RN"));
		m = suppl.next();
		assertEquals("78",m.getProp("_Name"));
		assertEquals("78",m.getProp("NSC"));
		assertEquals("6290-84-2",m.getProp("CAS_RN"));
	}
    @Test(expected=org.RDKit.GenericRDKitException.class)
	public void test1SDSupplierEOF() {
		File fileN = new File(baseTestPath, "NCI_aids.10.sdf");
		SDMolSupplier suppl = new SDMolSupplier(fileN.getPath());
		ROMol m;
		for (int i = 0; i < 10; i++)
			m = suppl.next();
		assertEquals( true,suppl.atEnd() );
		m = suppl.next();
	}

	@Test
	public void test2SmilesSupplier() {
		File fileN = new File(baseTestPath, "pgp_20.txt");
		SmilesMolSupplier suppl = 
			new SmilesMolSupplier(fileN.getPath(), "\t", 2, 1, true);

		ArrayList<ROMol> ms = new ArrayList<ROMol>();
		ROMol m;
		do {
			m = suppl.next();
			ms.add(m);
		}
		while (!suppl.atEnd());
		assertEquals( 20,ms.size() );
		suppl.reset();
		m = suppl.next();
		assertTrue (m.getProp("_Name").equals("ALDOSTERONE"));
		assertTrue (m.getProp("ID").equals("RD-PGP-0001"));
		m = suppl.next();
		assertTrue (m.getProp("_Name").equals("AMIODARONE"));
		assertTrue (m.getProp("ID").equals("RD-PGP-0002"));

		suppl.reset();
		assertEquals( 20 , suppl.length() );
	}

	@Test
	public void test3SmilesSupplier() throws Exception {
		String[] text = {
				"C1CC1,1",
				"CC(=O)O,3",
				"fail,4",
		"CCOC,5"};
		File f = File.createTempFile("smi", ".csv");
		String fileN = f.getPath();
		FileWriter fw = new FileWriter(f);
		for (String line : text)
			fw.write(line + "\n");
		fw.close();
		tmpFiles.add(fileN);
		SmilesMolSupplier suppl = new SmilesMolSupplier(fileN, ",", 0, 1, false);
		ArrayList<ROMol> ms = new ArrayList<ROMol>();
		ROMol m;
		do {
			m = suppl.next();
			if (m != null)
				ms.add(m);
		} while (!suppl.atEnd());
		assertEquals(3, ms.size());
	}


    //@Test
	public void test10Leak() {
            File fileN = new File(baseTestPath, "NCI_aids.10.sdf");
            SDMolSupplier suppl = new SDMolSupplier(fileN.getPath());
            Integer i;
            for(i=0;i<1000000;i++){
                ROMol m=suppl.next();
                m.delete();
                suppl.reset();
                if((i%1000)==0) System.err.printf("Done: %s\n",i);
            }
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.SuppliersTests");
	}

}
