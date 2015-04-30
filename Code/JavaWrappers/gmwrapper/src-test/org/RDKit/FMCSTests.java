/* 
 *
 *  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
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

public class FMCSTests extends GraphMolTest {
	private ROMol m;

	@Before public void setUp() {
	}

	@Test 
	public void test1Basics() {
		ROMol_Vect mols = new ROMol_Vect();
		mols.add(RWMol.MolFromSmiles("c1ccccc1OC"));
		mols.add(RWMol.MolFromSmiles("c1ccccc1C"));
		mols.add(RWMol.MolFromSmiles("c1c(C)cccc1OC"));
                
                MCSResult mcs=RDKFuncs.findMCS(mols);
                assertEquals(6,mcs.getNumAtoms());
                assertEquals(6,mcs.getNumBonds());
                assertEquals("[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",mcs.getSmartsString());
                assertEquals(false,mcs.getCanceled());

	}
	@Test 
	public void test2Basics() {
		ROMol_Vect mols = new ROMol_Vect();
		mols.add(RWMol.MolFromSmiles("c1ccccc1OC"));
		mols.add(RWMol.MolFromSmiles("c1ccccc1C"));
		mols.add(RWMol.MolFromSmiles("c1c(C)cccc1OC"));
		mols.add(RWMol.MolFromSmiles("C1CCCCC1C"));
                
                MCSResult mcs=RDKFuncs.findMCS(mols,true,1,60,false,false,false,false,false,
                                               AtomComparator.AtomCompareElements,
                                               BondComparator.BondCompareAny);
                assertEquals(6,mcs.getNumAtoms());
                assertEquals(6,mcs.getNumBonds());
                assertEquals("[#6]1:,-[#6]:,-[#6]:,-[#6]:,-[#6]:,-[#6]:,-1",mcs.getSmartsString());
                assertEquals(false,mcs.getCanceled());

	}
	@Test 
	public void test3Chirality() {
		ROMol_Vect mols = new ROMol_Vect();
        
		mols.add(RWMol.MolFromSmiles("C[C@H](F)CCl"));
		mols.add(RWMol.MolFromSmiles("C[C@H](F)C"));
                
                MCSResult mcs=RDKFuncs.findMCS_P(mols,"{\"MatchChiralTag\":true}");
                assertEquals(4,mcs.getNumAtoms());
                assertEquals(3,mcs.getNumBonds());
                //assertEquals("[#6]1:,-[#6]:,-[#6]:,-[#6]:,-[#6]:,-[#6]:,-1",mcs.getSmartsString());
                assertEquals(false,mcs.getCanceled());

	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.FMCSTests");
	}

}
