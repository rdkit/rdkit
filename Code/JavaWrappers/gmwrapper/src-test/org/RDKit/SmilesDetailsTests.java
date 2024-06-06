/* 
 * $Id: SmilesDetailsTests.java 134 2011-01-21 02:24:56Z bill.smith $
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

public class SmilesDetailsTests extends GraphMolTest {

	@Test
	public void testDetails() {
		ROMol mol;
		Atom a;
		String smi;

		// implicit/explicit H handling
		smi = "OC([OH])C[O-]";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(5, mol.getNumAtoms());
		a = mol.getAtomWithIdx(0);
		assertEquals(1, a.getImplicitValence());
		assertEquals(1, a.getExplicitValence());
		assertFalse(a.getNoImplicit());
		assertEquals(0, a.getFormalCharge());
		a = mol.getAtomWithIdx(2);
		assertEquals(0, a.getImplicitValence());
		assertEquals(2, a.getExplicitValence());
		assertTrue(a.getNoImplicit());
		assertEquals(0, a.getFormalCharge());
		a = mol.getAtomWithIdx(4);
		assertEquals(0, a.getImplicitValence());
		assertEquals(1, a.getExplicitValence());
		assertTrue(a.getNoImplicit());
		assertEquals(-1, a.getFormalCharge());

	}

	@Test
	public void testProblems() {
		ROMol mol;
		String smi;

		// ring closure handling with branches/fragments
		Int_Vect_Vect rings = new Int_Vect_Vect();
		smi = "C1(CC1CC1CC1)";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		int ringCount = RDKFuncs.findSSSR(mol, rings);
		assertEquals(2, ringCount);
		assertEquals(2, rings.size());
		assertEquals(3, rings.get(0).size());
		assertEquals(3, rings.get(1).size());

		// this is truly pathological, but both daylight
		// and chemdraw parse it properly
		smi = "C1.C1CC1CC1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		rings = new Int_Vect_Vect();
		ringCount = RDKFuncs.findSSSR(mol, rings);
		assertEquals(1, ringCount);
		assertEquals(1, rings.size());
		assertEquals(3, rings.get(0).size());

		// here's another stupid case that we need to handle:

		smi = "C1CC11CC1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		ringCount = RDKFuncs.findSSSR(mol, rings);
		assertEquals(2, ringCount);
		assertEquals(2, rings.size());
		assertEquals(3, rings.get(0).size());
		assertEquals(3, rings.get(1).size());

	}

	@Test
	public void testBasicCanon() {
		ROMol mol;
		String smi, refSmi;

		smi = "C1OCCCC1";
		mol = RWMol.MolFromSmiles(smi);
		refSmi = mol.MolToSmiles();

		smi = "C1COCCC1";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles();
		assertEquals(smi, refSmi);

		smi = "O1CCCCC1";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles();
		assertEquals(smi, refSmi);

		smi = "OC=CC";
		mol = RWMol.MolFromSmiles(smi);
		refSmi = mol.MolToSmiles();

		smi = "CC=CO";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles();
		assertEquals(smi, refSmi);

		smi = "C(C)=CO";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles();
		assertEquals(smi, refSmi);

		smi = "C(O)=CC";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles();
		assertEquals(smi, refSmi);

		// --- These are related to Issue 109

		smi = "C([H])Cl";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(2, mol.getNumAtoms());
		refSmi = mol.MolToSmiles();

		smi = "CCl";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles();
		assertEquals(smi, refSmi);

		// -- Issue 131
		smi = "P#[Ga]";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(2, mol.getNumAtoms());
		refSmi = mol.MolToSmiles();

		mol = RWMol.MolFromSmiles(refSmi);
		smi = mol.MolToSmiles();
		assertEquals(smi, refSmi);

		smi = "O=[Ba]";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(2, mol.getNumAtoms());
		refSmi = mol.MolToSmiles();

		mol = RWMol.MolFromSmiles(refSmi);
		smi = mol.MolToSmiles();
		assertEquals(smi, refSmi);

		// make sure empty molecules return empty SMILES:

		mol = new ROMol();
		smi = mol.MolToSmiles();
		assertEquals("", smi);

	}

	@Test
	public void testStereochem() {
		ROMol mol;
		String smi, refSmi, cip;

		smi = "F[C@](Cl)(Br)I";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(0).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);
		refSmi = mol.MolToSmiles(true);

		smi = "F[C@](Br)(I)Cl";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(0).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "F[C@](I)(Cl)Br";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "Cl[C@](Br)(F)I";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "Cl[C@](F)(I)Br";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "I[C@](F)(Br)Cl";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "I[C@](Br)(Cl)F";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "F[C@@](Br)(Cl)I";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "F[C@@](Cl)(I)Br";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "Cl[C@@](Br)(I)F";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "Cl[C@@](F)(Br)I";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "[C@@](Cl)(F)(Br)I";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(0).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(0).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(0).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "F[C@H](Cl)Br";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("R", cip);
		refSmi = mol.MolToSmiles(true);

		smi = "Br[C@H](F)Cl";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		smi = mol.MolToSmiles(true);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("R", cip);
		assertEquals(refSmi, smi);

		smi = "Br[C@]([H])(F)Cl";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("R", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "Br[C@](F)(Cl)[H]";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("R", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "Br[C@]1(F)(Cl).[H]1";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("R", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "Br[C@H]1Cl.F1";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("R", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "Br[C@]12Cl.F2.[H]1";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("R", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "Br[C@]21Cl.F1.[H]2";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("R", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "[C@@H](Br)(F)Cl";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(0).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(0).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(0).getProp("_CIPCode");
		assertEquals("R", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		smi = "[H][C@@](Br)(F)Cl";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(0).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(0).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(0).getProp("_CIPCode");
		assertEquals("R", cip);
		smi = mol.MolToSmiles(true);
		assertEquals(refSmi, smi);

		// an additional set of test cases from the Chirality notes document.
		// one can never have too many tests of this stuff.

		smi = "F[C@]([H])(O)C";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);

		smi = "F[C@]1([H])OC1";
		mol = RWMol.MolFromSmiles(smi);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);

		smi = "F[C@H](O)C";
		mol = RWMol.MolFromSmiles(smi);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);

		smi = "F[C@@H]1OC1";
		mol = RWMol.MolFromSmiles(smi);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("S", cip);

		smi = "[C@](F)([H])(O)C";
		mol = RWMol.MolFromSmiles(smi);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(0).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(0).getProp("_CIPCode");
		assertEquals("S", cip);

		smi = "[C@@]1(F)([H])OC1";
		mol = RWMol.MolFromSmiles(smi);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(0).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(0).getProp("_CIPCode");
		assertEquals("S", cip);

		smi = "[C@@H](F)(O)C";
		mol = RWMol.MolFromSmiles(smi);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(0).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(0).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals("C[C@@H](O)F", smi);
		smi = mol.MolToSmiles(true, false, 0);
		assertEquals("[C@H](C)(O)F", smi);

		smi = "[C@@H]1(F)OC1";
		mol = RWMol.MolFromSmiles(smi);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(0).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(0).getProp("_CIPCode");
		assertEquals("S", cip);
		smi = mol.MolToSmiles(true);
		assertEquals("F[C@H]1CO1", smi);
		smi = mol.MolToSmiles(true, false, 0);
		assertEquals("[C@H]1(F)CO1", smi);

		smi = "C1O[C@H]1F";
		mol = RWMol.MolFromSmiles(smi);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(2).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(2).getProp("_CIPCode");
		assertEquals("S", cip);

		smi = "C1O[C@@]1([H])F";
		mol = RWMol.MolFromSmiles(smi);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(2).hasProp("_CIPCode"));
		cip = mol.getAtomWithIdx(2).getProp("_CIPCode");
		assertEquals("S", cip);

		// -----------------------------------
		// test some double-bond containing molecules:

		// -- cis --

		smi = "F\\C=C/Br";
		mol = RWMol.MolFromSmiles(smi);
		refSmi = mol.MolToSmiles(true);

		mol = RWMol.MolFromSmiles(refSmi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "Br\\C=C/F";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "Br/C=C\\F";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "F/C=C\\Br";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		// -- trans --

		smi = "F\\C=C\\Br";
		mol = RWMol.MolFromSmiles(smi);
		refSmi = mol.MolToSmiles(true);

		mol = RWMol.MolFromSmiles(refSmi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "Br\\C=C\\F";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "Br/C=C/F";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "F/C=C/Br";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		// -- more complex --

		smi = "F\\C=C(/Cl)\\Br";
		mol = RWMol.MolFromSmiles(smi);
		refSmi = mol.MolToSmiles(true);

		mol = RWMol.MolFromSmiles(refSmi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "F/C=C(\\Cl)/Br";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "F/C=C(\\Cl)Br";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "F/C=C(Cl)/Br";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		// -- combine chirality with cis/trans --

		smi = "F[C@H](Cl)\\C=C(/F)";
		mol = RWMol.MolFromSmiles(smi);
		refSmi = mol.MolToSmiles(true);

		mol = RWMol.MolFromSmiles(refSmi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "F[C@H](Cl)/C=C(\\F)";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "Cl[C@@H](F)/C=C(\\F)";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "Cl[C@@H](F)\\C=C(/F)";
		mol = RWMol.MolFromSmiles(smi);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

	}

	@Test
	public void testIssue127() {
		ROMol mol, mol2;
		String smi, refSmi, tempStr;

		smi = "Cl[C@]12[Si]C(C2)O1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);

		// first roundtrip the non-chiral SMILES:
		refSmi = mol.MolToSmiles();
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		tempStr = mol2.MolToSmiles();
		assertEquals(tempStr, refSmi);

		// now do the true SMILES:
		refSmi = mol.MolToSmiles(true);
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		tempStr = mol2.MolToSmiles(true);
		assertEquals(tempStr, refSmi);

	}

	@Test
	public void testIssue143() {
		ROMol mol;
		String smi, refSmi;

		smi = "C[C@](C)(C)C";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		refSmi = mol.MolToSmiles(true);
		assertEquals("CC(C)(C)C", refSmi);

		smi = "CC[C@](C)(C)C=O";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		refSmi = mol.MolToSmiles(true);
		assertEquals("CCC(C)(C)C=O", refSmi);

	}

	@Test
	public void testIssue151() {
		ROMol mol, mol2;
		String smi, refSmi;

		smi = "C1S[C@H]1O";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(0).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(1).getChiralTag());
		assertTrue(mol.getAtomWithIdx(2).getChiralTag() != Atom.ChiralType.CHI_UNSPECIFIED);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(2).getChiralTag());

		refSmi = mol.MolToSmiles(true);
		assertEquals("O[C@H]1CS1", refSmi);
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		smi = mol2.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "F[C@@H]1O[C@H](Cl)S1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(0).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(2).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(4).getChiralTag());

		assertTrue(mol.getAtomWithIdx(1).getChiralTag() != Atom.ChiralType.CHI_UNSPECIFIED);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		assertTrue(mol.getAtomWithIdx(3).getChiralTag() != Atom.ChiralType.CHI_UNSPECIFIED);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(3).getChiralTag());

		refSmi = mol.MolToSmiles(true);
		assertEquals("F[C@@H]1O[C@H](Cl)S1", refSmi);
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		smi = mol2.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "Cl[C@@H]1S[C@@H](O1)F";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(0).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(2).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(4).getChiralTag());

		assertTrue(mol.getAtomWithIdx(1).getChiralTag() != Atom.ChiralType.CHI_UNSPECIFIED);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		assertTrue(mol.getAtomWithIdx(3).getChiralTag() != Atom.ChiralType.CHI_UNSPECIFIED);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(3).getChiralTag());

		refSmi = mol.MolToSmiles(true);
		assertEquals("F[C@@H]1O[C@H](Cl)S1", refSmi);
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		smi = mol2.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "Cl[C@@H]1O[C@H](F)S1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(0).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(2).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(4).getChiralTag());

		assertTrue(mol.getAtomWithIdx(1).getChiralTag() != Atom.ChiralType.CHI_UNSPECIFIED);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		assertTrue(mol.getAtomWithIdx(3).getChiralTag() != Atom.ChiralType.CHI_UNSPECIFIED);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(3).getChiralTag());

		refSmi = mol.MolToSmiles(true);
		assertEquals("F[C@H]1O[C@@H](Cl)S1", refSmi);
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		smi = mol2.MolToSmiles(true);
		assertEquals(smi, refSmi);

	}

	@Test
	public void testIssue153() {
		String code;
		ROMol mol, mol2;
		String smi, refSmi;

		smi = "C1(O[C@H]12)S2";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(0).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(1).getChiralTag());
		assertTrue(mol.getAtomWithIdx(2).getChiralTag() != Atom.ChiralType.CHI_UNSPECIFIED);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(2).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(2).hasProp("_CIPCode"));
		code = mol.getAtomWithIdx(2).getProp("_CIPCode");
		assertEquals("S", code);

		refSmi = mol.MolToSmiles(true);
		assertEquals("O1C2S[C@H]12", refSmi);
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		smi = mol2.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "C1(O[C@H]21)S2";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(0).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(1).getChiralTag());
		assertTrue(mol.getAtomWithIdx(2).getChiralTag() != Atom.ChiralType.CHI_UNSPECIFIED);
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(2).getChiralTag());
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(2).hasProp("_CIPCode"));
		code = mol.getAtomWithIdx(2).getProp("_CIPCode");
		assertEquals("R", code);

		refSmi = mol.MolToSmiles(true);
		assertEquals("O1C2S[C@@H]12", refSmi);
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		smi = mol2.MolToSmiles(true);
		assertEquals(smi, refSmi);

	}

	@Test
	public void testIssue157() {
		String code;
		ROMol mol, mol2;
		String smi, refSmi;
		smi = "O[C@](C)(Cl)[C@@](O)(Cl)C";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Atom.ChiralType.CHI_UNSPECIFIED, mol.getAtomWithIdx(0).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CCW, mol.getAtomWithIdx(1).getChiralTag());
		assertEquals(Atom.ChiralType.CHI_TETRAHEDRAL_CW, mol.getAtomWithIdx(4).getChiralTag());

		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		code = mol.getAtomWithIdx(1).getProp("_CIPCode");
		assertEquals("R", code);
		assertTrue(mol.getAtomWithIdx(1).hasProp("_CIPCode"));
		code = mol.getAtomWithIdx(4).getProp("_CIPCode");
		assertEquals("S", code);

		refSmi = mol.MolToSmiles(true);
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		smi = mol2.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "Cl[C@@](C)1CC[C@@](C)(C1)Cl";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		refSmi = mol.MolToSmiles(true);
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		smi = mol2.MolToSmiles(true);
		assertEquals(smi, refSmi);

		smi = "[H][C@@]12CC(CO1)CN2";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		RDKFuncs.assignStereochemistry(mol);
		smi = mol.getAtomWithIdx(0).getProp("_CIPCode");
		assertEquals("S", smi);
		refSmi = mol.MolToSmiles(true);

		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		smi = mol2.MolToSmiles(true);

		assertEquals(smi, refSmi);

		smi = "[H][C@@]12C[14C@@](C=C1)(C3C2C(NC3=O)=O)[H]";
		// smi="C1=C[C@@H]2C[C@H]1C1C(=O)NC(=O)C21";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		RDKFuncs.assignStereochemistry(mol);
		smi = mol.getAtomWithIdx(0).getProp("_CIPCode");
		assertEquals("R", smi);
		smi = mol.getAtomWithIdx(2).getProp("_CIPCode");
		assertEquals("S", smi);
		refSmi = mol.MolToSmiles(true);
		mol2 = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol2);
		smi = mol2.MolToSmiles(true);

		assertEquals(smi, refSmi);

	}

	@Test
	public void testIssue159() {
		ROMol mol;
		String smi, refSmi;

		smi = "C/C=C/O";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);

		assertEquals(Bond.BondStereo.STEREONONE, mol.getBondWithIdx(0).getStereo());
		assertEquals(Bond.BondStereo.STEREOE, mol.getBondWithIdx(1).getStereo());
		refSmi = mol.MolToSmiles(true);

		smi = "C(\\C)=C/O";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);

		assertEquals(Bond.BondStereo.STEREONONE, mol.getBondWithIdx(0).getStereo());
		assertEquals(Bond.BondStereo.STEREOE, mol.getBondWithIdx(1).getStereo());
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);
	}

	@Test
	public void testIssue175() {
		ROMol mol;
		String smi;

		smi = "Cl\\C=C1.F/1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Bond.BondStereo.STEREOE, mol.getBondWithIdx(1).getStereo());

		smi = "Cl\\C=C1CN/1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Bond.BondStereo.STEREOE, mol.getBondWithIdx(1).getStereo());
	}

	@Test
	public void testIssue176() {
		ROMol mol;
		String smi;
		smi = "C1CC1C1CC1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(7, mol.getNumBonds());

		smi = "C1CC1C1CC-1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(7, mol.getNumBonds());

		smi = "C1CC1C1CC=1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(7, mol.getNumBonds());

		smi = "C1CC1C=1CC1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(7, mol.getNumBonds());
	}

	// Keep this -- first test of BondIterator
	@Test
	public void testIssue184() {
		ROMol mol;
		String smi, refSmi;

		smi = "C1NC(Cl)C(=N\\O)/C1=N\\O";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Bond.BondType.DOUBLE, mol.getBondWithIdx(4).getBondType());
		assertEquals(Bond.BondStereo.STEREOZ, mol.getBondWithIdx(4).getStereo());
		assertEquals(Bond.BondType.DOUBLE, mol.getBondWithIdx(7).getBondType());
		assertEquals(Bond.BondStereo.STEREOZ, mol.getBondWithIdx(7).getStereo());
		refSmi = mol.MolToSmiles(true);

		mol = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);

		for (BondIterator bondIt = mol.beginBonds(); bondIt.ne(mol.endBonds()); bondIt.next()) {
			Bond b = bondIt.getBond();
			if (b.getBondType() == Bond.BondType.DOUBLE) {
				assertEquals(Bond.BondStereo.STEREOZ, b.getStereo());
			}
		}

	}

	@Test
	public void testIssue185() {
		ROMol mol;
		String smi, refSmi;

		// start with a simple E/Z handling case with branches:
		smi = "C(/C)=N/O";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Bond.BondType.DOUBLE, mol.getBondWithIdx(1).getBondType());
		assertEquals(Bond.BondStereo.STEREOZ, mol.getBondWithIdx(1).getStereo());
		refSmi = mol.MolToSmiles(true, false, 0); // (1,0,0);

		assertEquals("C(\\C)=N\\O", refSmi);

		// make sure we can round-trip:
		mol = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol);
		assertEquals(Bond.BondType.DOUBLE, mol.getBondWithIdx(1).getBondType());
		assertEquals(Bond.BondStereo.STEREOZ, mol.getBondWithIdx(1).getStereo());

	}

	@Test
	public void testIssue191() {
		ROMol mol;
		String smi, refSmi;
		int numE = 0;

		smi = "C2=NNC(N=C2)=N\\N=C\\c1ccccc1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(Bond.BondType.DOUBLE, mol.getBondWithIdx(7).getBondType());
		assertEquals(Bond.BondStereo.STEREOE, mol.getBondWithIdx(7).getStereo());
		refSmi = mol.MolToSmiles(true);

		mol = RWMol.MolFromSmiles(refSmi);
		assertNotNull(mol);
		numE = 0;
		for (BondIterator bondIt = mol.beginBonds(); bondIt.ne(mol.endBonds()); bondIt.next()) {
			Bond b = bondIt.getBond();
			if (b.getBondType() == Bond.BondType.DOUBLE) {
				assertTrue(b.getStereo() != Bond.BondStereo.STEREOZ);
				if (b.getStereo() == Bond.BondStereo.STEREOE) {
					numE++;
				}
			}
		}
		assertEquals(1, numE);
		smi = mol.MolToSmiles(true);
		assertEquals(smi, refSmi);
	}

	@Test
	public void testIssue256() {
		ROMol mol;
		Bond bond;
		String smi;


		smi = "C1CC[C+]1=1CCC1";
		mol = RWMol.MolFromSmiles(smi, 0, false);
		assertNotNull(mol);
		bond = mol.getBondBetweenAtoms(3, 0);
		assertNotNull(bond);
		assertEquals(Bond.BondType.SINGLE, bond.getBondType());
		bond = mol.getBondBetweenAtoms(3, 6);
		assertNotNull(bond);
		assertEquals(Bond.BondType.DOUBLE, bond.getBondType());

		smi = "C1CC[C+]=11CCC1";
		mol = RWMol.MolFromSmiles(smi, 0, false);
		assertNotNull(mol);
		bond = mol.getBondBetweenAtoms(3, 0);
		assertNotNull(bond);
		assertEquals(Bond.BondType.DOUBLE, bond.getBondType());
		bond = mol.getBondBetweenAtoms(3, 6);
		assertNotNull(bond);
		assertEquals(Bond.BondType.SINGLE, bond.getBondType());

	}

	@Test
	public void testIssue266() {
		RWMol mol;
		String smi;

		smi = "c1ccccc1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles();
		assertEquals("c1ccccc1", smi);

		RDKFuncs.Kekulize(mol);
		smi = mol.MolToSmiles();
		assertEquals("C1=CC=CC=C1", smi);

		smi = "c1ccccc1c1ccccc1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles();
		assertEquals("c1ccc(-c2ccccc2)cc1", smi);

		RDKFuncs.Kekulize(mol);
		smi = mol.MolToSmiles();
		assertEquals("C1=CC=C(C2=CC=CC=C2)C=C1", smi);

	}

	@Test
	public void testRootedAt() {
		RWMol mol;
		String smi;

		smi = "CN(C)C";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(false, false, -1);
		assertEquals("CN(C)C", smi);
		smi = mol.MolToSmiles(false, false, 1);
		assertEquals("N(C)(C)C", smi);
		smi = mol.MolToSmiles(false, false, 2);
		assertEquals("CN(C)C", smi);

	}

	@Test
	public void testIsotopes() {
		RWMol mol;
		String smi;

		smi = "C[13C](C)(C)C";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		assertEquals(0.001, mol.getAtomWithIdx(1).getMass(), 13.0034);
                assertEquals(mol.getAtomWithIdx(1).getIsotope(), 13);
		smi = mol.MolToSmiles(false);
		assertEquals("CC(C)(C)C", smi);
		smi = mol.MolToSmiles(true);
		assertEquals("C[13C](C)(C)C", smi);

	}

	@Test
	public void testBug1670149() {
		RWMol mol;
		String smi;

		smi = "C1[NH2+]CCC1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(false, false, -1);
		assertEquals("C1CC[NH2+]C1", smi);

		mol.getAtomWithIdx(1).setNumExplicitHs(0);
		mol.getAtomWithIdx(1).setNoImplicit(false);
		mol.getAtomWithIdx(1).updatePropertyCache();
		assertEquals(2, mol.getAtomWithIdx(1).getNumImplicitHs());
		smi = mol.MolToSmiles(false, false, -1);
		assertEquals("C1CC[NH2+]C1", smi);

	}

	@Test
	public void testBug1719046() {
		RWMol mol;
		String smi;

		smi = "Cl[CH]1CCCCC1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(false, false, -1);
		assertEquals("ClC1CCCCC1", smi);

		smi = "Cl[C@H]1CCCCC1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(false, false, -1);
		assertEquals("ClC1CCCCC1", smi);

		smi = "Cl[C@H]1C(Br)CCCC1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(false, false, -1);
		assertEquals("ClC1CCCCC1Br", smi);

		smi = "[CH]1=[CH][CH]=[CH][CH]=[CH]1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(false, false, -1);
		assertEquals("c1ccccc1", smi);

		smi = "c1ccccn1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(false, false, -1);
		assertEquals("c1ccncc1", smi);

		smi = "C1=CNC=C1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(false, false, -1);
		assertEquals("c1cc[nH]c1", smi);

		smi = "[CH]1=[CH][NH][CH]=[CH]1";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(false, false, -1);
		assertEquals("c1cc[nH]c1", smi);

	}

	@Test
	public void testBug1842174() {

		RWMol mol;
		String smi;

		smi = "F/C=N/Cl";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(true, false, -1);

		assertEquals("F/C=N/Cl", smi);

		smi = mol.MolToSmiles(true, false, 1);

		assertEquals("C(\\F)=N/Cl", smi);

		smi = "C(\\C=C\\F)=C(/Cl)Br";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi = mol.MolToSmiles(true, false, -1);

		assertEquals("F/C=C/C=C(/Cl)Br", smi);

	}

	@Test
	public void testBug1844617() {

		RWMol mol;
		String smi, smi2;
		String label;

		smi = "O=C1C2OCC[C@@]22C(CC1)CNCC2";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(6).hasProp("_CIPCode"));
		label = mol.getAtomWithIdx(6).getProp("_CIPCode");
		assertEquals("S", label);

		smi = mol.MolToSmiles(true);

		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi2 = mol.MolToSmiles(true);

		assertEquals(smi2, smi);

		smi = "O=C1CC[C@@]2(O)[C@@H]3N(C)CC[C@]22[C@H]1OC[C@H]2CC3";
		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		RDKFuncs.assignStereochemistry(mol);
		assertTrue(mol.getAtomWithIdx(4).hasProp("_CIPCode"));
		label = mol.getAtomWithIdx(4).getProp("_CIPCode");
		assertEquals("S", label);
		assertTrue(mol.getAtomWithIdx(6).hasProp("_CIPCode"));
		label = mol.getAtomWithIdx(6).getProp("_CIPCode");
		assertEquals("R", label);
		assertTrue(mol.getAtomWithIdx(11).hasProp("_CIPCode"));
		label = mol.getAtomWithIdx(11).getProp("_CIPCode");
		assertEquals("S", label);
		assertTrue(mol.getAtomWithIdx(12).hasProp("_CIPCode"));
		label = mol.getAtomWithIdx(12).getProp("_CIPCode");
		assertEquals("R", label);
		assertTrue(mol.getAtomWithIdx(15).hasProp("_CIPCode"));
		label = mol.getAtomWithIdx(15).getProp("_CIPCode");
		assertEquals("S", label);

		smi = mol.MolToSmiles(true);

		mol = RWMol.MolFromSmiles(smi);
		assertNotNull(mol);
		smi2 = mol.MolToSmiles(true);

		assertEquals(smi2, smi);
	}

	@Test
	public void testRingStereochem() {

		RWMol m;
		String smi;

		smi = "C[C@H]1CC[C@@H](C)CC1";
		m = RWMol.MolFromSmiles(smi);
		assertNotNull(m);
		assertEquals(8, m.getNumAtoms());

		smi = m.MolToSmiles(true);
		//assertTrue(m.hasProp("_ringStereoWarning"));

		smi = m.MolToSmiles(false);
		//assertTrue((!m.hasProp("_ringStereoWarning")));

	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.SmilesDetailsTests");
	}

}
