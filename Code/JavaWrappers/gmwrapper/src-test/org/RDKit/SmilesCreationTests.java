/* 
 * $Id: SmilesCreationTests.java 131 2011-01-20 22:01:29Z ebakke $
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

import org.junit.Test;

public class SmilesCreationTests extends GraphMolTest {

	static String[] goodSMILES =
	{
		// We probably don't need to test so many SMILES
		// in the Java wrapper since we're relying on the C++
		// code to be working, but it doesn't hurt.
		"C1CC2C1CC2",
		"c1cccn(=O)c1",
		"C",
		"CC",
		"C-C",
		"C=C",
		"[CH2+]C[CH+2]",
		"C1CC1",
		"C1CC=1",
		"C=1CC1",
		"C=C-O",
		"C1CC1",
		"C1NC1",
		"C1=CC1",
		"C1CCC1",
		"CC(C)CC",
		"CC(=O)O",
		"C1C(=O)C1",
		"C1C(N)C1",
		"CC(O)C",
		"OC=CCC",
		"CC([O-])O",
		"C1CC2C1CC2",
		"Cl/C=C/Cl",
		"Cl/C=C\\Cl",
		"Cl/C=C/Cl",
		"Cl/C=C\\Cl",
		"C1CC.CC1",
		"C1C(C2CC2).C2CC2C1",
		"[Na+].[Cl-].[NH4+].[Cl-]",
		"C[35Cl]",
		"C%10CC%10",
		"[H][H]",
		"[H+]",
		"C[N+](=O)[O-]",
		"N1C(=N)SC=C1",
		"[O-][N+](=O)C1=CNC(=N)S1",
		"CN(=O)=O",
		"C1=CC=C[N+]([O-])=C1",
		"C1=CC=CN(=O)=C1",		
		"CC(=CO)C",
		"CCC",
		"C1CC1"
	};
	
	static String[] whitespaceSMILES =
	{
		// test whitespace tolerance:
		"  C1=CC=CN(=O)=C1",
		"C1=CC=CN(=O)=C1  ",
		"  C1=CC=CN(=O)=C1  ",
		"\tC1=CC=CN(=O)=C1\r\n"
	};
		
	static String[] dummyAtomSMILES =
	{
		// test dummy atoms:
		"c1ccccc1[*]",
		"c1ccccc1[1*]",
		"S1cccc1",
		"*1ccccc1",
		"C1=CC=CC=C1",
		"*1=CC=CC=C1",
		"*1*cccc1",
		"*1**ccc1"
	};
	
	static String[] aromaticSeAndTeSmiles =
	{
		// test aromatic se and te:
		"c1ccc[se]1",
		"c1ccc[te]1"
	};
	
	@Test
	public void testGoodSMILES() {
		for (String smi : goodSMILES) {
			testGoodSMILES(smi);
		}
	}

	@Test
	public void testWhitespaceSMILES() {
		for (String smi : whitespaceSMILES) {
			testGoodSMILES(smi);
		}
	}

	@Test
	public void testWDummyAtomSMILES() {
		for (String smi : dummyAtomSMILES) {
			testGoodSMILES(smi);
		}
	}

	@Test
	public void testAromaticSeAndTeSMILES() {
		for (String smi : aromaticSeAndTeSmiles) {
			testGoodSMILES(smi);
		}
	}

	@Test
	public void testZerosAsRingIndices() {
		// test zeros as ring indices, issue 2690982:
		testGoodSMILES("C0CC0");
	}

	@Test
	public void testCanonization() {
		// test canonization error, issue 3018558:
		testGoodSMILES("C/C(/C=C2\\Sc1ccc(cc1N\\2C))=C5\\SC4=NccN4C\\5=O");
	}
	
	
	/*
	 * Test creation failures
	 */
	@Test
	public void testBadSMILESBadBranchSyntax() {
		assertNull( RWMol.MolFromSmiles("CC=(CO)C") );
	}
	@Test
	public void testBadSMILESBadRingSytax() {
		assertNull( RWMol.MolFromSmiles("C1CC"));
	}
	@Test(expected=MolSanitizeException.class)
	public void testBadSMILESBadAromaticity() {
		RWMol.MolFromSmiles("Ccc");
	}
	@Test
	public void testBadSMILESNonsenseInput() {
		assertNull( RWMol.MolFromSmiles("fff"));
	}
	@Test(expected=MolSanitizeException.class)
	public void testBadSMILESBadValence() {
		RWMol.MolFromSmiles("N(=O)(=O)=O");
	}
	@Test
	public void testBadSMILESIncorrectHandlingOfZeroAsRingClosureDigit_1() {
		// part of sf.net issue 2525792
		assertNull( RWMol.MolFromSmiles("C=0"));
	}
	@Test
	public void testBadSMILESIncorrectHandlingOfZeroAsRingClosureDigit_2() {
		// part of sf.net issue 2525792
		assertNull( RWMol.MolFromSmiles("C0"));
	}
	@Test
	public void testBadSMILESIncorrectHandlingOfZeroAsRingClosureDigit_3() {
		// part of sf.net issue 2525792
		assertNull( RWMol.MolFromSmiles("C-0"));
	}
	@Test
	public void testBadSMILESIncorrectHandlingOfZeroAsRingClosureDigit_4() {
		// part of sf.net issue 2525792
		assertNull( RWMol.MolFromSmiles("C+0"));
	}

	private void testGoodSMILES(String smi) {
		ROMol mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);

		long nAts = mol.getNumAtoms();
		assertFalse(nAts == 0);
		String smi2 = ((RWMol) mol).MolToSmiles();
		ROMol mol2 = RWMol.MolFromSmiles(smi2);
		if (mol2 != null) {
			assertEquals(nAts, mol2.getNumAtoms());
		}
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.SmilesCreationTests");
	}
}
