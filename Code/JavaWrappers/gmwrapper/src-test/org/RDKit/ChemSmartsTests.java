/* 
 * $Id: ChemSmartsTests.java 131 2011-01-20 22:01:29Z ebakke $
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

import org.junit.*;

public class ChemSmartsTests extends GraphMolTest {
	private ROMol m;

	@Before
	public void setUp() {
		String rdpath = System.getenv("RDBASE");
		if (rdpath == null)
			org.junit.Assert.fail("No definition for RDBASE");
		File base = new File(rdpath);
		File testFile = new File(base, "rdkit" + File.separator + "Chem"
				+ File.separator + "test_data" + File.separator + "quinone.mol");
		String fn = testFile.getAbsolutePath();
		m = RWMol.MolFromMolFile(fn);
	}

	// testing molecule
	@Test
	public void testMol() {
		assertEquals("bad nAtoms", 8, m.getNumAtoms());
	}

	// testing smarts match
	@Test
	public void testMatch() {
		ROMol p = RWMol.MolFromSmarts("CC(=O)C");
		Match_Vect_Vect matches = m.getSubstructMatches(p);
		assertEquals("bad match count:  " + matches.size(), 2, matches.size());
		for (int i = 0; i < matches.size(); i++)
			assertEquals("bad match (" + i + ")", 4, matches.get(i).size());
	}

	// test atom order in smarts match
	@Test
	public void testOrder() {
		ROMol p = RWMol.MolFromSmarts("CC(=[O,N])C");
		Match_Vect_Vect matches = m.getSubstructMatches(p);
		Match_Vect match = matches.get(0);
		String atoms = "";
		for (int i = 0; i < match.size(); i++)
			atoms += m.getAtomWithIdx(match.get(i).getSecond()).getSymbol();
		assertEquals("bad atom ordering:  " + atoms, "CCOC", atoms);

	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.ChemSmartsTests");
	}

}
