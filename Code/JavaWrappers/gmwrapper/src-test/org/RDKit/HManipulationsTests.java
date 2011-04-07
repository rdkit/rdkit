/* 
 * $Id: HManipulationsTests.java 131 2011-01-20 22:01:29Z ebakke $
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

public class HManipulationsTests extends GraphMolTest {

	@Before
	public void setUp() {
	}
	@Test
	public void testExplicitHs() {
		ROMol m = RWMol.MolFromSmiles("CC(=O)[OH]");
		assertEquals(m.getNumAtoms(), 4);

		ROMol m2 = m.addHs(true);
		assertEquals(m2.getNumAtoms(), 5);
		m2 = m2.removeHs(true);
		assertEquals(m2.getNumAtoms(), 5);
	}
	@Test
	public void testImplicitHs() {
		ROMol m = RWMol.MolFromSmiles("CC(=O)[OH]");
		ROMol m2 = m.addHs(true);
		m2 = m2.removeHs(true);

		m2 = m2.removeHs(false);
		assertEquals(m2.getNumAtoms(), 4);

		m2 = m.addHs(false);
		assertEquals(m2.getNumAtoms(), 8);
	}
	@Test
	public void testExplicitAndImplicitHs() {
		ROMol m = RWMol.MolFromSmiles("CC(=O)[OH]");
		ROMol m2 = m.addHs(true);
		assertEquals(m2.getNumAtoms(), 5);
		m2 = m2.removeHs(true);
		assertEquals(m2.getNumAtoms(), 5);
		m2 = m2.removeHs(false);
		assertEquals(m2.getNumAtoms(), 4);
	}
	@Test
	public void testMergingQueryHs() {
		// Create a molecule without sanitizing -- expect
		// that the hydrogen is kept even though it should be
		// implicit.
		ROMol m = RWMol.MolFromSmiles("CC[H]", 0, false);
		assertEquals(m.getNumAtoms(), 3);
		ROMol m2 = m.mergeQueryHs();
		assertEquals(m2.getNumAtoms(), 2);
		assertTrue(m2.getAtomWithIdx(1).hasQuery());
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.HManipulationsTests");
	}

}
