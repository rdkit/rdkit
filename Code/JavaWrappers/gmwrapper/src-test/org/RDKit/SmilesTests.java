/* 
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

public class SmilesTests extends GraphMolTest {

	private void testSpellings(String smi, String[] spellings) {
		ROMol m = RWMol.MolFromSmiles(smi);
		String canSmi = m.MolToSmiles();
		for (String spelling : spellings) {
			m = RWMol.MolFromSmiles(spelling);
			assertNotNull("Can't parse " + spelling, m);

			String trySmi = m.MolToSmiles();
			assertEquals("Non-canonical:  mol " + spelling + " gave " + trySmi + "(should be "
					+ canSmi + ")", canSmi, trySmi);

			m = RWMol.MolFromSmiles(trySmi);
			String trySmi2 = m.MolToSmiles();
			assertEquals("Non-canonical:  mol " + spelling + " gave " + trySmi2 + "(should be "
					+ canSmi + ") on second pass", canSmi, trySmi2);
		}
	}

	// testing first batch of linear mols
	@Test
	public void testLinear1() {
		testSpellings("O=CCO", new String[] { "OCC=O", "C(O)C=O", "C(C=O)O", "C(CO)=O" });

		testSpellings("OCC(C=C)CCC(C#N)CC", new String[] { "C=CC(CO)CCC(C#N)CC",
				"C(CO)(C=C)CCC(CC)C#N", "C(CO)(C=C)CCC(C#N)CC", "C(C=C)(CO)CCC(C#N)CC",
				"C(C=C)(CO)CCC(CC)C#N" });

		testSpellings("[Se]=CCO", new String[] { "OCC=[Se]", "C(O)C=[Se]", "C(C=[Se])O","C(CO)=[Se]" });
	}

	// testing first batch of rings
	@Test
	public void testRings1() {
		testSpellings("C1OCCCC1", new String[] { "O1CCCCC1", "C1COCCC1", "C1CCOCC1", "C1CCCOC1","C1CCCCO1"});
		testSpellings("CC1=CCCCC1", new String[] { "C1=C(C)CCCC1", "C1CC=C(C)CC1"});
		testSpellings("CC1C=CCCC1", new String[] { "C1=CC(C)CCC1", "C1CC=CC(C)C1"});
	}

	// testing second batch of rings
	@Test
	public void testRings2() {
		testSpellings("c1c(cc2nc3cc(ccc3cc2c1))", new String[] { "c1ccc2cc3ccccc3nc2c1",
				"c1ccc2nc3ccccc3cc2c1", "c1c2nc3ccccc3cc2ccc1"});
		testSpellings("Cc1ccc2nc3ccccc3cc2c1", new String[] { "c1ccc2nc3ccc(C)cc3cc2c1"});
		testSpellings("c1c(C)cc2nc3ccccc3cc2c1", new String[] { "c1ccc2nc3cc(C)ccc3cc2c1"});
	}

	// testing molecules which have been problematic
	@Test
	public void testProblems() {
		testSpellings("[Al+2]CCC",
				new String[] { "CCC[Al+2]", "C(C)(C[Al+2])"});
		testSpellings("C(=O)(Cl)CC(=O)Cl", 
				new String[] { "ClC(CC(Cl)=O)=O", "C(Cl)(=O)CC(=O)Cl","C(Cl)(=O)CC(Cl)=O"});
		testSpellings("C(=O)(Cl)c1ccc(C(=O)Cl)cc1",
				new String[] { "O=C(Cl)c1ccc(cc1)C(Cl)=O","C(Cl)(=O)C1=CC=C(C=C1)C(Cl)=O", "ClC(=O)c1ccc(cc1)C(=O)Cl"});
		testSpellings("[N+](=O)([O-])c1ccc([N+](=O)[O-])cc1",
				new String[] { "[N+]([O-])(=O)C1=CC=C(C=C1)[N+](=O)[O-]","O=[N+1]([O-1])c1ccc(cc1)[N+1]([O-1])=O", "[O-1][N+1](=O)c1ccc(cc1)[N+1]([O-1])=O"});
		testSpellings("Oc1c3c(cc(c1)S(=O)(=O)O)cc(NC(=O)c2ccccc2)cc3",
				new String[] { "C1=C(C2=C(C=C1S(O)(=O)=O)C=C(C=C2)NC(C3=CC=CC=C3)=O)O", "O=S(=O)(O)c1cc(O)c2ccc(NC(=O)c3ccccc3)cc2c1", "OS(=O)(=O)c1cc(O)c2ccc(NC(=O)c3ccccc3)cc2c1"});
		testSpellings("C",
				new String[] { "C"});
		testSpellings("C(Cl)(Br)(F)CC(Cl)(Br)(F)",
				new String[] { "C(Cl)(F)(Br)CC(F)(Cl)(Br)","C(Cl)(Br)(F)CC(Cl)(F)(Br)", "C(F)(Br)(Cl)CC(Br)(Cl)(F)","C(C(Cl)(Br)(F))C(F)(Cl)Br"});
	}

	// testing tricky (high-symmetry) molecules
	@Test
	public void testHighSymmetry() {
		testSpellings("CC(C)CC", new String[] { "CCC(C)C"});
		testSpellings("C1CCCC1CCC", new String[] { "CCCC1CCCC1"});
		testSpellings("C1(C)CC(C)CCC1", new String[] { "CC1CCCC(C)C1"});
		testSpellings("CCC1CCCCC1CC", new String[] { "CCC1CCCCC1CC"});
		testSpellings("CCC1CC(CC)CCC1", new String[] { "CCC1CCCC(CC)C1"});
		testSpellings("C1CCCCC1CC(CC)CC", new String[] { "CCC(CC)CC1CCCCC1"});
		testSpellings("C1CCCC2C1CC(CC)CC2", new String[] { "CCC1CCC2CCCCC2C1"});
		testSpellings("CC1CCCC2C1C(C)CCC2", new String[] { "CC1CCCC2CCCC(C)C12"});
		testSpellings("C2CCC1CCC(C)C12", new String[] { "CC1CCC2CCCC12"});
		testSpellings("CC(C)CCCC(C)C", new String[] { "CC(CCCC(C)C)C"});
	}

	// EXPECT FAILURES -> testing molecules which are known to fail
	@Test
	public void testFailures() {
		testSpellings("C13C6C1C2C4C2C3C5C4C56", 
				new String[] { "C45C1C6C3C6C5C4C2C3C12","C45C2C6C3C6C5C4C1C3C12"});
	}

	@Test
	public void testReplacements() {
            String_String_Map repls=new String_String_Map();
            repls.set("{X}","OC1CC1");
            RWMol nmol = RWMol.MolFromSmiles("c1ccccc1{X}",0,true,repls); 
            String nsmi = RDKFuncs.MolToSmiles(nmol, true);
            String expected="c1ccc(OC2CC2)cc1";
            assertEquals("bad smiles: "+nsmi+"!="+expected,nsmi,expected);
	}

        @Test
	public void testRankAtoms(){
	    //Need a molecule to canonicalise
	    // expected ordering here: [11, 8, 3, 5, 0, 9, 7, 10, 6, 1, 4, 2]
	    ROMol m1 = RWMol.MolFromSmiles("C(CO)(C=C)CCC(CC)C#N");

	    // same molecule, different atom ordering:
	    // expected ordering here: [11, 5, 0, 8, 3, 9, 7, 10, 4, 2, 6, 1] 
	    ROMol m2 = RWMol.MolFromSmiles("C(C=C)(CO)CCC(C#N)CC");

	    UInt_Vect ranks1 = new UInt_Vect();
	    m1.rankMolAtoms(ranks1);
	    assertEquals("Wrong size ranks - " + ranks1.size() + " != " + 
	  		m1.getNumAtoms(), ranks1.size(), m1.getNumAtoms());

	    UInt_Vect ranks2 = new UInt_Vect();
	    m2.rankMolAtoms(ranks2);
	    assertEquals("Wrong size ranks - " + ranks2.size() + " != " + 
	  		m2.getNumAtoms(), ranks2.size(), m2.getNumAtoms());

	    Match_Vect_Vect matches = m1.getSubstructMatches(m2);
	    assertEquals("bad matches size: "+matches.size(),matches.size(),1);
	    Match_Vect match = matches.get(0);
	    assertEquals("bad match size: "+match.size(),match.size(),m1.getNumAtoms());
	    for(int i=0;i<match.size();i++){
	 	assertEquals("bad rank: "+match.get(i)+" "+ranks1+" "+ranks2,
			ranks1.get(match.get(i).getSecond()),
			ranks2.get(match.get(i).getFirst()));
	    }
	    
	    m1.delete();
	    m2.delete();
	    ranks1.delete();
	    ranks2.delete();
	    matches.delete();
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.SmilesTests");
	}

}
