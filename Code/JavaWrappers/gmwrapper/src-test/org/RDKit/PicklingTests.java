/* 
 * $Id: PicklingTests.java 131 2011-01-20 22:01:29Z ebakke $
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
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import org.junit.*;

public class PicklingTests extends GraphMolTest {

	private ROMol mol1;

	@Before
	public void setUp() {
		String smiles = "c1ccccc1";
		mol1 = RWMol.MolFromSmiles(smiles);
	}

	@Test
	public void testBasicPickling() {
		Int_Vect pkl = mol1.ToBinary();
		ROMol m1 = ROMol.MolFromBinary(pkl);
		assertEquals(6, m1.getNumAtoms());
		assertEquals(6, m1.getNumBonds());
	}

	@Test
	public void testBasicReactionPickling() {
		ChemicalReaction rxn;
		rxn = ChemicalReaction.ReactionFromSmarts("[OH][C:1]=[O:2].[N!H0:3]>>[N:3][C:1]=[O:2]");
		assertEquals(2, rxn.getNumReactantTemplates());
		assertEquals(1, rxn.getNumProductTemplates());
		Int_Vect pkl = rxn.ToBinary();
		ChemicalReaction rxn2 = ChemicalReaction.RxnFromBinary(pkl);
		assertEquals(2, rxn2.getNumReactantTemplates());
		assertEquals(1, rxn2.getNumProductTemplates());
	}

	@Test
	public void testIssue164() {
		String smi = "NCCCCC(NC(C(C(C)C)NC(=O)C1CCCN1C(C(N)CCCNC(=N)N)=O)=O)C(NC(C(C)C)C(NC(Cc2ccc(O)cc2)C(=O)N6C(C(NC(CC(N)=O)C(NCC(NC(C)C(NC(CCC(O)=O)C(NC(CC(O)=O)C(NC(C(NC(CO)C(NC(C)C(NC(CCC(O)=O)C(NC(C)C(NC(Cc3ccccc3)C(=O)N5C(C(NC(CC(C)C)C(NC(C(NC(C(O)=O)Cc4ccccc4)=O)CCC(O)=O)=O)=O)CCC5)=O)=O)=O)=O)=O)CCC(O)=O)=O)=O)=O)=O)=O)=O)CCC6)=O)=O";
		ROMol m1 = RWMol.MolFromSmiles(smi);
		Int_Vect pickle;
		ROMol m2;
		pickle = m1.ToBinary();
		m2 = ROMol.MolFromBinary(pickle);

		assertEquals(m2.getNumAtoms(), m1.getNumAtoms());

		// the issue had to do with number of atoms, so let's make an enormous
		// molecule and try again:
		RWMol m3 = RWMol.MolFromSmiles(smi);
		m3.insertMol(m1);
		m3.insertMol(m1);
		m3.sanitizeMol();

		pickle = m3.ToBinary();
		ROMol m4 = ROMol.MolFromBinary(pickle);

		assertEquals(m4.getNumAtoms(), m3.getNumAtoms());
		assertEquals(m4.getNumBonds(), m3.getNumBonds());

	}


	@Test
	public void testIssue219() {

		String smi = "CC";
		ROMol m1 = RWMol.MolFromSmiles(smi);

		Int_Vect pickle = m1.ToBinary();
		ROMol m2 = ROMol.MolFromBinary(pickle);
		assertEquals(m2.getNumAtoms(),m1.getNumAtoms());

		Conformer conf = new Conformer(2);
		conf.setId(23);
		m1.addConformer(conf);
		pickle = m1.ToBinary();
		m2 = ROMol.MolFromBinary(pickle);
		assertEquals(m2.getNumAtoms(),m1.getNumAtoms());
		assertEquals(23,m2.getConformer().getId());

	}

	@Test
	public void testIssue220() {

		String smi = "N/N=C/C";
		ROMol m1 = RWMol.MolFromSmiles(smi);
		assertEquals(m1.getBondWithIdx(1).getStereo(),Bond.BondStereo.STEREOE);
		Int_Vect pickle = m1.ToBinary();
		ROMol m2 = ROMol.MolFromBinary(pickle);
		assertEquals(m2.getNumAtoms(), m1.getNumAtoms());
		assertEquals(m2.getNumBonds(), m1.getNumBonds());
		assertEquals(m2.getBondWithIdx(1).getStereo(), Bond.BondStereo.STEREOE);
		for (int i = 0; i < m1.getNumBonds(); ++i) {
			assertEquals(m1.getBondWithIdx(i).getStereo(), m2.getBondWithIdx(i).getStereo());
			Int_Vect s1 = m1.getBondWithIdx(i).getStereoAtoms();
			Int_Vect s2 = m2.getBondWithIdx(i).getStereoAtoms();
			assertEquals(s2.size(), s1.size());
			for (int j = 0; j < s1.size(); j++)
				assertEquals(s2.get(j), s1.get(j));
		}
	}

	@Test
	public void testToFromByteArray() throws IOException {
		String smi = "CN(C)c1ccc2c(=O)cc[nH]c2c1";
		String pklFileName = "quinolone.pkl";
		{
			ROMol mol = RWMol.MolFromSmiles(smi);
			byte[] pkl = mol.toByteArray();
			FileOutputStream pklOutStream = null;
			try {
				pklOutStream = new FileOutputStream(pklFileName);
				pklOutStream.write(pkl);
			} finally {
				if (pklOutStream != null) {
					pklOutStream.close();
				}
			}
			mol.delete();
		}
		{
			FileInputStream pklInStream = null;
			byte[] pkl = null;
			File pklInFile = new File(pklFileName);
			try {
				pklInStream = new FileInputStream(pklInFile);
				pkl = new byte[(int)pklInFile.length()];
				assertEquals(pklInStream.read(pkl), pkl.length);
			} finally {
				if (pklInStream != null) {
					pklInStream.close();
				}
			}
			ROMol mol = ROMol.fromByteArray(pkl);
			assertEquals(mol.MolToSmiles(), smi);
			mol.delete();
		}
	}

	@Test
	public void testPickleProperties() {
		ROMol mol = RWMol.MolFromSmiles("c1ccccc1[C@](F)(Cl)Br");
		mol.setProp("foo", "bar");
		mol.setProp("_MolFileChiralFlag", "1");
		{ 
			Int_Vect pkl = mol.ToBinary();
			ROMol mol2 = ROMol.MolFromBinary(pkl);
			assertFalse(mol2.hasProp("_MolFileChiralFlag"));
			assertFalse(mol2.hasProp("foo"));
		}
		{ 
			Int_Vect pkl = mol.ToBinary(PropertyPickleOptions.AllProps.swigValue());
			ROMol mol2 = ROMol.MolFromBinary(pkl);
			assertTrue(mol2.hasProp("_MolFileChiralFlag"));
			assertTrue(mol2.hasProp("foo"));
		}
		{ 
			long val = RDKFuncs.getDefaultPickleProperties();
			RDKFuncs.setDefaultPickleProperties(PropertyPickleOptions.AllProps.swigValue());
			Int_Vect pkl = mol.ToBinary();
			RDKFuncs.setDefaultPickleProperties(val);
			ROMol mol2 = ROMol.MolFromBinary(pkl);
			assertTrue(mol2.hasProp("_MolFileChiralFlag"));
			assertTrue(mol2.hasProp("foo"));
		}

	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.PicklingTests");
	}
}
