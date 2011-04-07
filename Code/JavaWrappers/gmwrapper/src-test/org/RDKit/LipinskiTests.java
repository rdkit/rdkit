/* 
 * $Id: LipinskiTests.java 131 2011-01-20 22:01:29Z ebakke $
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

public class LipinskiTests extends GraphMolTest {
	private File dataFile;	

	@Before
	public void setUp() {
		File base = getRdBase();
		dataFile = new File(base, "Data" + File.separator + "NCI" + File.separator + 
		"first_200.props.sdf");
	}
	@Test
	public void test1() {
		SDMolSupplier suppl = new SDMolSupplier(dataFile.getPath());
		int idx = 1;
		while (!suppl.atEnd()) {
			long calc;
			long orig;

			ROMol m = suppl.next();
			calc = RDKFuncs.calcNumHBD(m);
			orig = Integer.parseInt(m.getProp("NUM_HDONORS"));
			assertTrue("bad num h donors for mol " + idx + " ("
					+ m.getProp("SMILES") + "): " + calc + " != " + orig,
					calc==orig); 

			calc = RDKFuncs.calcNumHBA(m);
			orig = Integer.parseInt(m.getProp("NUM_HACCEPTORS"));
			assertTrue("bad num h acceptors for mol " + idx + " ("
					+ m.getProp("SMILES") + "): " + calc + " != " + orig,
					calc==orig); 

			calc = RDKFuncs.calcNumHeteroatoms(m);
			orig = Integer.parseInt(m.getProp("NUM_HETEROATOMS"));
			assertTrue("bad num heteroatoms for mol " + idx + " ("
					+ m.getProp("SMILES") + "): " + calc + " != " + orig,
					calc==orig); 

			calc = RDKFuncs.calcNumRotatableBonds(m);
			orig = Integer.parseInt(m.getProp("NUM_ROTATABLEBONDS"));
			assertTrue("bad num rotatable bonds for mol " + idx + " ("
					+ m.getProp("SMILES") + "): " + calc + " != " + orig,
					calc==orig); 
			idx += 1;
		}
	}

	// testing a problem with acceptor definition
	@Test
	public void testIssue2183420() {
		assertEquals(1,RDKFuncs.calcNumHBA(RWMol.MolFromSmiles("NC")));
		assertEquals(1,RDKFuncs.calcNumHBA(RWMol.MolFromSmiles("CNC")));
		assertEquals(1,RDKFuncs.calcNumHBA(RWMol.MolFromSmiles("CN(C)C")));
		assertEquals(1,RDKFuncs.calcNumHBA(RWMol.MolFromSmiles("NC(=O)")));
		assertEquals(1,RDKFuncs.calcNumHBA(RWMol.MolFromSmiles("NC(=O)C")));
		assertEquals(1,RDKFuncs.calcNumHBA(RWMol.MolFromSmiles("CNC(=O)")));
		assertEquals(1,RDKFuncs.calcNumHBA(RWMol.MolFromSmiles("CNC(=O)C")));
		assertEquals(2,RDKFuncs.calcNumHBA(RWMol.MolFromSmiles("O=CNC(=O)C")));
		assertEquals(2,RDKFuncs.calcNumHBA(RWMol.MolFromSmiles("O=C(C)NC(=O)C")));
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.LipinskiTests");
	}

}
