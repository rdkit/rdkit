/* 
 * $Id$
 *
 *  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
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

public class ChemTransformsTests extends GraphMolTest {
	private ROMol m;

	@Before public void setUp() {
	}

	@Test 
	public void test1BRICS() {
		ROMol mol = RWMol.MolFromSmiles("c1ccccc1OC");
                ROMol nmol=RDKFuncs.fragmentOnBRICSBonds(mol);
		assertEquals( 8,mol.getNumAtoms());
		assertEquals( 10,nmol.getNumAtoms());
                assertEquals("[3*]OC.[16*]c1ccccc1",nmol.MolToSmiles(true));

	}
	@Test 
	public void test2BRICS() {
		ROMol mol = RWMol.MolFromSmiles("OC(C)=CC");
                ROMol nmol=RDKFuncs.fragmentOnBRICSBonds(mol);
		assertEquals( 5,mol.getNumAtoms());
		assertEquals( 7,nmol.getNumAtoms());
                assertEquals("[7*]=CC.[7*]=C(C)O",nmol.MolToSmiles(true));

	}
	@Test 
	public void test3BRICS() {
		ROMol mol = RWMol.MolFromSmiles("CCCOCCC(=O)c1ccccc1");
                ROMol nmol=RDKFuncs.fragmentOnBRICSBonds(mol);
		assertEquals( 14,mol.getNumAtoms());
		assertEquals( 20,nmol.getNumAtoms());
                assertEquals("[3*]O[3*].[4*]CCC.[4*]CCC([6*])=O.[16*]c1ccccc1",nmol.MolToSmiles(true));

	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.ChemTransformsTests");
	}

}
