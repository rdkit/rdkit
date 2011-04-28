/* 
* $Id: ChemAtomTests.java 131 2011-01-20 22:01:29Z ebakke $
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

import org.junit.*;

public class ChemAtomTests extends GraphMolTest {

	private ROMol m;

	@Before
	public void setUp() {
		String smiles = "CC(=O)CCSC";
		m = RWMol.MolFromSmiles(smiles);
	}

	// A new method introduced in the ROMol.i wrapper
	@Test
	public void testGetAtoms() {
		Atom_Vect atoms = m.getAtoms();
		assertEquals(7L, atoms.size());

	}

	@Test
	public void test1Implicit() {
		// Test implicit valence
		Atom a = m.getAtoms().get(0);
		int iV = a.getImplicitValence();
		assertEquals(3, iV);
		assertEquals(0, m.getAtomWithIdx(1).getImplicitValence());
		assertEquals(0, m.getAtomWithIdx(2).getImplicitValence());
		assertEquals(2, m.getAtomWithIdx(3).getImplicitValence());
	}

	// Test getBonds
	@Test
	public void testGetBonds() {
		Bond_Vect bonds = m.getAtomWithIdx(1).getBonds();
		assertEquals(3, bonds.size());
		assertEquals(1, m.getAtomWithIdx(2).getBonds().size());
	}

	// Testing iteration over an atom's bonds
	@Test
	public void test2BondIter() {
		Atom a = m.getAtomWithIdx(1);
		Bond_Vect bonds = a.getBonds();
		assertEquals(3, bonds.size());
		for (int bond_num = 0; bond_num < bonds.size(); bond_num++) {
			assertNotNull(bonds.get(bond_num));
		}
	}

	// testing GetBondBetweenAtoms(idx,idx)
	@Test
	public void test3GetBond() {
		Bond b = m.getBondBetweenAtoms(1, 2);
		assertEquals("GetBondBetweenAtoms failed", Bond.BondType.DOUBLE, b.getBondType());
	}

	// Testing atomic props
	@Test
	public void test4Props() {
		Atom a = m.getAtomWithIdx(1);
		assertEquals("C",a.getSymbol());
		assertEquals(6, a.getAtomicNum());
		assertEquals(0, a.getFormalCharge());
		assertEquals(3, a.getDegree());
		assertEquals(0, a.getImplicitValence());
		assertEquals(4, a.getExplicitValence());
	}

	// Testing setting atomic props and generic exception handler
	@Test
	public void test5Setters() {
		Atom a = new Atom(6);
		assertEquals("C",a.getSymbol());
		assertEquals( 6,a.getAtomicNum() );
		a.setFormalCharge(1);
		assertEquals( 1,a.getFormalCharge() );
	}
	
	@Test(expected=GenericRDKitException.class)
	public void testExceptionHandler() {
		Atom a = new Atom(6);	
		a.getImplicitValence();
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.ChemAtomTests");
	}

}
