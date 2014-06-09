/* 
 * $Id: ChemBondTests.java 131 2011-01-20 22:01:29Z ebakke $
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
public class ChemBondTests extends GraphMolTest {
	private ROMol m;

	@Before public void setUp() {
		String smiles="CCCC1=CC=C1";
		m = RWMol.MolFromSmiles(smiles);
	}

	// testing GetBond
	@Test 
	public void test1Get() {
		assertNotNull( m.getBondBetweenAtoms(0, 1) );
	}

	// Testing setting bond props
	@Test
	public void test2Setters() {
		Bond b = m.getBondBetweenAtoms(0,1);
		assertEquals(Bond.BondType.SINGLE, b.getBondType());
		b.setBondDir(Bond.BondDir.BEGINWEDGE);
		assertEquals(Bond.BondDir.BEGINWEDGE,m.getBondBetweenAtoms(0,1).getBondDir());
	}

	// Testing bond props
	@Test
	public void test3Props() {
		Bond b = m.getBondBetweenAtoms(0,1);
		assertEquals(Bond.BondType.SINGLE,b.getBondType());
		assertEquals(m.getAtomWithIdx(0).getIdx(),b.getBeginAtom().getIdx());
		assertEquals(0,b.getBeginAtomIdx());
		assertEquals(m.getAtomWithIdx(1).getIdx(),b.getEndAtom().getIdx());
		assertEquals(1,b.getEndAtomIdx());
	}

	// Testing more bond props
	@Test
	public void test4Props2() {
		Bond b = m.getBondBetweenAtoms(3,4);
		assertEquals(Bond.BondType.DOUBLE,b.getBondType());
		Bond b2 = m.getBondBetweenAtoms(1,2);
		assertEquals(Bond.BondType.SINGLE,b2.getBondType());
		assertTrue (b.getIsConjugated());
		assertFalse (b2.getIsConjugated());
	}

	// Testing more bond props
	@Test
	public void test5BondIters() {
	    // Traversing through bond list of a molecule
	    BondIterator it = m.beginBonds();
	    int iCount = 0;
 
	    while ( it.ne(m.endBonds()) ) {
		final Bond bond = it.getBond();
		final Bond.BondType type = bond.getBondType();
		final Atom atomStart = bond.getBeginAtom();
		final Atom atomEnd = bond.getEndAtom();
		// some silly tests
		assertEquals(atomStart.getIdx(),bond.getBeginAtomIdx());
		assertEquals(atomEnd.getIdx(),bond.getEndAtomIdx());

		it = it.next();
		iCount++;
	    }
	    assertEquals(iCount,m.getNumBonds());
	}

    public static void main(String args[]) {
	org.junit.runner.JUnitCore.main("org.RDKit.ChemBondTests");
    }

}
