/* 
 * $Id: ChemTests.java 131 2011-01-20 22:01:29Z ebakke $
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

public class AssignChiralityTests extends GraphMolTest {
	
	@Test
	public void testAssignStereochemistry() {
		
		//Atom 1 is possibly chiral, atom 5 is definitely chiral
		ROMol m = RWMol.MolFromSmiles("CC(Cl)(F)C[C@H](O)C");
		
		//Try the old method...
		RDKFuncs.assignStereochemistry(m, true, true, true);
		Atom at = m.getAtomWithIdx(1);
		assertTrue("Atom 1 is STEREOANY", at.hasProp("_ChiralityPossible") 
			   && at.getChiralTag() == Atom.ChiralType.CHI_UNSPECIFIED);
		at = m.getAtomWithIdx(5);
		assertTrue("Atom 5 is Has Sterero", at.getChiralTag() == Atom.ChiralType.CHI_TETRAHEDRAL_CCW 
			   || at.getChiralTag() == Atom.ChiralType.CHI_TETRAHEDRAL_CW);
	
		//Now try new method with default for optional boolean
		Int_Vect centres = new Int_Vect();
		RDKFuncs.assignStereochemistry(m, true, true, true, centres);
		at=m.getAtomWithIdx(1);
		assertTrue("Atom 1 is STEREOANY", at.hasProp("_ChiralityPossible") 
			   && at.getChiralTag() == Atom.ChiralType.CHI_UNSPECIFIED);
		at=m.getAtomWithIdx(5);
		assertTrue("Atom 5 is Has Sterero", at.getChiralTag() == Atom.ChiralType.CHI_TETRAHEDRAL_CCW 
			   || at.getChiralTag() == Atom.ChiralType.CHI_TETRAHEDRAL_CW);
		assertTrue("Atoms 1 and 5 are listed", centres.get(0) == 1 && centres.get(1) == 5);
		centres.delete();
		
		//Now try new method with onlyStereoAny=false
		centres = new Int_Vect();
		RDKFuncs.assignStereochemistry(m, true, true, true, centres, false);
		at = m.getAtomWithIdx(1);
		assertTrue("Atom 1 is STEREOANY", at.hasProp("_ChiralityPossible") 
			   && at.getChiralTag() == Atom.ChiralType.CHI_UNSPECIFIED);
		at=m.getAtomWithIdx(5);
		assertTrue("Atom 5 is Has Sterero", at.getChiralTag() == Atom.ChiralType.CHI_TETRAHEDRAL_CCW 
			   || at.getChiralTag() == Atom.ChiralType.CHI_TETRAHEDRAL_CW);
		assertTrue("Atoms 1 and 5 are listed", centres.get(0) == 1 && centres.get(1) == 5);
		centres.delete();
		
		//Now try new method with onlyStereoAny=false
		centres = new Int_Vect();
		RDKFuncs.assignStereochemistry(m, true, true, true, centres, true);
		at = m.getAtomWithIdx(1);
		assertTrue("Atom 1 is STEREOANY", at.hasProp("_ChiralityPossible") && 
			   at.getChiralTag() == Atom.ChiralType.CHI_UNSPECIFIED);
		at = m.getAtomWithIdx(5);
		assertTrue("Atom 5 is Has Sterero", at.getChiralTag() == Atom.ChiralType.CHI_TETRAHEDRAL_CCW 
			   || at.getChiralTag() == Atom.ChiralType.CHI_TETRAHEDRAL_CW);
		assertTrue("Atoms 1 is listed", centres.get(0)==1 && centres.size()==1);
		centres.delete();
		
		m.delete();
	}

	
	@Test
	public void testFindPotentialStereobonds(){
		ROMol m = RWMol.MolFromSmiles("CC=CC");
		
		//Old method...
		RDKFuncs.findPotentialStereoBonds(mol, true);
		assertTrue("Bond 1 is possible stereo", 
			   m.getBondWithIdx(1).getStereo() == Bond.BondStereo.STEREONONE);
		
		//New method...
		Int_Vect stereoBonds = new Int_Vect();
		RDKFuncs.findPotentialStereoBonds(mol, true, stereoBonds);
		assertTrue("Bond 1 is possible stereo", 
			   m.getBondWithIdx(1).getStereo() == Bond.BondStereo.STEREONONE 
			   && stereoBonds.get(0) == 1 && stereoBonds.size() == 1);
		stereoBonds.delete();
		m.delete();
	}
	
	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.AssignChiralityTests");
	}

}
