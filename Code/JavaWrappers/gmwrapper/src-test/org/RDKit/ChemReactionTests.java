/* 
 * $Id: ChemReactionTests.java 131 2011-01-20 22:01:29Z ebakke $
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

import java.io.*;
import java.util.ArrayList;

import org.junit.*;

public class ChemReactionTests extends GraphMolTest {
	private ArrayList<String> tmpFiles = new ArrayList<String>();
	private File baseTestPath;
	@Before 
	public void setUp() {
		File base = getRdBase();
		baseTestPath = new File(base, "Code" + File.separator + "GraphMol" + File.separator + 
				"ChemReactions" + File.separator + "testData");
	}

	@After
	public void tearDown() throws Exception {
		for (String fn : tmpFiles) {
			new File(fn).delete();
		}
	}

	@Test
	public void test1Basics() {
		ChemicalReaction rxn = new ChemicalReaction();
		assertEquals( 0,rxn.getNumReactantTemplates() );
		assertEquals( 0,rxn.getNumProductTemplates() );

		ROMol r1 = RWMol.MolFromSmarts("[C:1](=[O:2])O");
		rxn.addReactantTemplate(r1);
		assertEquals( 1,rxn.getNumReactantTemplates() );

		r1 = RWMol.MolFromSmarts("[N;!$(N-C=O):3]");
		rxn.addReactantTemplate(r1);
		assertEquals( 2,rxn.getNumReactantTemplates() );

		r1 = RWMol.MolFromSmarts("[C:1](=[O:2])[N:3]");
		rxn.addProductTemplate(r1);
		assertEquals( 1,rxn.getNumProductTemplates() );

		ROMol_Vect reacts = new ROMol_Vect();
		for (String react : new String[] {"C(=O)O","N"})
			reacts.add(RWMol.MolFromSmiles(react));

		// Need this init call since the reaction has been built
		// "manually"
		rxn.initReactantMatchers();
		ROMol_Vect_Vect prods = rxn.runReactants(reacts);;
		assertEquals( 1,prods.size() );
		assertEquals( 1,prods.get(0).size() );
		assertEquals( 3,prods.get(0).get(0).getNumAtoms() );
                assertEquals( true, RDKFuncs.isMoleculeReactantOfReaction(rxn,reacts.get(0)) );
                assertEquals( false, RDKFuncs.isMoleculeReactantOfReaction(rxn,prods.get(0).get(0)) );
                assertEquals( true, RDKFuncs.isMoleculeProductOfReaction(rxn,prods.get(0).get(0)) );
                assertEquals( false, RDKFuncs.isMoleculeProductOfReaction(rxn,reacts.get(0)) );
	}

	@Test
	public void test2DaylightParser() {
		ChemicalReaction rxn = 
			ChemicalReaction.ReactionFromSmarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]");
		assertNotNull(rxn );;
		assertEquals( 2,rxn.getNumReactantTemplates() );
		assertEquals( 1,rxn.getNumProductTemplates() );
		assertTrue(rxn.getImplicitPropertiesFlag());

		ROMol_Vect reacts = new ROMol_Vect();
		for (String react : new String[] {"C(=O)O","N"})
			reacts.add(RWMol.MolFromSmiles(react));		
		// Do not need the initReactantMatchers call here
		ROMol_Vect_Vect prods = rxn.runReactants(reacts);;
		assertEquals( 1,prods.size() );
		assertEquals( 1,prods.get(0).size() );
		assertEquals( 3,prods.get(0).get(0).getNumAtoms() );

		reacts.clear();
		for (String react : new String[] {"CC(=O)OC","CN"})
			reacts.add(RWMol.MolFromSmiles(react));		
		// Do not need the initReactantMatchers call here
		prods = rxn.runReactants(reacts);;
		assertEquals( 1,prods.size() );
		assertEquals( 1,prods.get(0).size() );
		assertEquals( 5,prods.get(0).get(0).getNumAtoms() );
	}
	@Test
	public void test3MDLParsers() throws Exception {
		File rxnFile = new File(baseTestPath, "AmideBond.rxn");
		String fileN = rxnFile.getPath();
		ChemicalReaction rxn = ChemicalReaction.ReactionFromRxnFile(fileN);
		assertNotNull(rxn );;
		assertEquals( 2,rxn.getNumReactantTemplates() );
		assertEquals( 1,rxn.getNumProductTemplates() );
		assertFalse(rxn.getImplicitPropertiesFlag());

		ROMol_Vect reacts = new ROMol_Vect();
		for (String react : new String[] {"C(=O)O","N"})
			reacts.add(RWMol.MolFromSmiles(react));		
		// Do not need the initReactantMatchers call here
		ROMol_Vect_Vect prods = rxn.runReactants(reacts);;
		assertEquals( 1,prods.size() );
		assertEquals( 1,prods.get(0).size() );
		assertEquals( 3,prods.get(0).get(0).getNumAtoms() );

		String rxnBlock = "";
		BufferedReader f = new BufferedReader(new FileReader(rxnFile));
		String line;
		try{
			while ((line = f.readLine()) != null)
				rxnBlock += line + "\n"; // ok only for small numbers of lines
		} finally {
			f.close();
		}
		rxn = ChemicalReaction.ReactionFromRxnBlock(rxnBlock);
		assertNotNull(rxn );;
		assertEquals( 2,rxn.getNumReactantTemplates() );
		assertEquals( 1,rxn.getNumProductTemplates() );

		reacts.clear();
		for (String react : new String[] {"C(=O)O","N"})
			reacts.add(RWMol.MolFromSmiles(react));		
		// Do not need the initReactantMatchers call here
		prods = rxn.runReactants(reacts);;
		assertEquals( 1,prods.size() );
		assertEquals( 1,prods.get(0).size() );
		assertEquals( 3,prods.get(0).get(0).getNumAtoms() );
	}

	@Test (expected=ChemicalReactionParserException.class)
	public void test4ErrorHandling_1() {
		ChemicalReaction.ReactionFromSmarts("[C:1](=[O:2])Q.[N:3]>>[C:1](=[O:2])[N:3]");
	}
	@Test (expected=ChemicalReactionParserException.class)
	public void test4ErrorHandling_2() {
		ChemicalReaction.ReactionFromSmarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]Q");
	}
	@Test (expected=ChemicalReactionParserException.class)
	public void test4ErrorHandling_3() {
		ChemicalReaction.ReactionFromSmarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]>>CC");
	}
	@Test (expected=ChemicalReactionParserException.class)
	public void test4ErrorHandling_4() {
		String rxnBlock = 
			"$RXN" + "\n" +  
			"" + "\n" +  
			"      ISIS     082120061354" + "\n" +  
			"" + "\n" +  
			"  3  1" + "\n" +  
			"$MOL" + "\n" +  
			"" + "\n" +  
			"  -ISIS-  08210613542D" + "\n" +  
			"" + "\n" +  
			"  3  2  0  0  0  0  0  0  0  0999 V2000" + "\n" +  
			"   -1.4340   -0.6042    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0" + "\n" +  
			"   -0.8639   -0.9333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0" + "\n" +  
			"   -1.4340    0.0542    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0" + "\n" +  
			"  1  2  1  0  0  0  0" + "\n" +  
			"  1  3  2  0  0  0  0" + "\n" +  
			"M  END" + "\n" +  
			"$MOL" + "\n" +  
			"" + "\n" +  
			"  -ISIS-  08210613542D" + "\n" +  
			"" + "\n" +  
			"  1  0  0  0  0  0  0  0  0  0999 V2000" + "\n" +  
			"    2.2125   -0.7833    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0" + "\n" +  
			"M  END" + "\n" +  
			"$MOL" + "\n" +  
			"" + "\n" +  
			"  -ISIS-  08210613542D" + "\n" +  
			"" + "\n" +  
			"  3  2  0  0  0  0  0  0  0  0999 V2000" + "\n" +  
			"    9.5282   -0.8083    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0" + "\n" +  
			"    8.9579   -0.4792    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0" + "\n" +  
			"    8.9579    0.1792    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0" + "\n" +  
			"  1  2  1  0  0  0  0" + "\n" +  
			"  2  3  2  0  0  0  0" + "\n" +  
			"M  END";
		ChemicalReaction.ReactionFromRxnBlock(rxnBlock);
	}
	@Test (expected=ChemicalReactionParserException.class)
	public void test4ErrorHandling_5() {
		String rxnBlock = 
			"$RXN" + "\n" + 
			"" + "\n" + 
			"      ISIS     082120061354" + "\n" + 
			"" + "\n" + 
			"  2  1" + "\n" + 
			"$MOL" + "\n" + 
			"" + "\n" + 
			"  -ISIS-  08210613542D" + "\n" + 
			"" + "\n" + 
			"  4  2  0  0  0  0  0  0  0  0999 V2000" + "\n" + 
			"   -1.4340   -0.6042    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0" + "\n" + 
			"   -0.8639   -0.9333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0" + "\n" + 
			"   -1.4340    0.0542    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0" + "\n" + 
			"  1  2  1  0  0  0  0" + "\n" + 
			"  1  3  2  0  0  0  0" + "\n" + 
			"M  END" + "\n" + 
			"$MOL" + "\n" + 
			"" + "\n" + 
			"  -ISIS-  08210613542D" + "\n" + 
			"" + "\n" + 
			"  1  0  0  0  0  0  0  0  0  0999 V2000" + "\n" + 
			"    2.2125   -0.7833    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0" + "\n" + 
			"M  END" + "\n" + 
			"$MOL" + "\n" + 
			"" + "\n" + 
			"  -ISIS-  08210613542D" + "\n" + 
			"" + "\n" + 
			"  3  2  0  0  0  0  0  0  0  0999 V2000" + "\n" + 
			"    9.5282   -0.8083    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0" + "\n" + 
			"    8.9579   -0.4792    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0" + "\n" + 
			"    8.9579    0.1792    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0" + "\n" + 
			"  1  2  1  0  0  0  0" + "\n" + 
			"  2  3  2  0  0  0  0" + "\n" + 
			"M  END" + "\n" ;
		ChemicalReaction.ReactionFromRxnBlock(rxnBlock);
	}
	@Test (expected=ChemicalReactionParserException.class)
	public void test4ErrorHandling_6() {
		String rxnBlock = 	
			"$RXN" + "\n" + 
			"" + "\n" + 
			"      ISIS     082120061354" + "\n" + 
			"" + "\n" + 
			"  2  1" + "\n" + 
			"$MOL" + "\n" + 
			"" + "\n" + 
			"  -ISIS-  08210613542D" + "\n" + 
			"" + "\n" + 
			"  3  2  0  0  0  0  0  0  0  0999 V2000" + "\n" + 
			"   -1.4340   -0.6042    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0" + "\n" + 
			"   -0.8639   -0.9333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0" + "\n" + 
			"   -1.4340    0.0542    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0" + "\n" + 
			"  1  2  1  0  0  0  0" + "\n" + 
			"  1  3  2  0  0  0  0" + "\n" + 
			"M  END" + "\n" + 
			"$MOL" + "\n" + 
			"" + "\n" + 
			"  -ISIS-  08210613542D" + "\n" + 
			"" + "\n" + 
			"  1  0  0  0  0  0  0  0  0  0999 V2000" + "\n" + 
			"    2.2125   -0.7833    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0" + "\n" + 
			"M  END" + "\n" + 
			"$MOL" + "\n" + 
			"" + "\n" + 
			"  -ISIS-  08210613542D" + "\n" + 
			"" + "\n" + 
			"  3  1  0  0  0  0  0  0  0  0999 V2000" + "\n" + 
			"    9.5282   -0.8083    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0" + "\n" + 
			"    8.9579   -0.4792    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0" + "\n" + 
			"    8.9579    0.1792    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0" + "\n" + 
			"  1  2  1  0  0  0  0" + "\n" + 
			"  2  3  2  0  0  0  0" + "\n" + 
			"M  END" + "\n";
		ChemicalReaction.ReactionFromRxnBlock(rxnBlock);
	}

	@Test
	public void test5Validation() {
		ChemicalReaction rxn = ChemicalReaction.ReactionFromSmarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]");
		assertNotNull(rxn);
		Int_Pair validationResults = rxn.validateReaction();
		assertEquals(0, validationResults.getFirst());
		assertEquals(0, validationResults.getSecond());

		rxn = ChemicalReaction.ReactionFromSmarts("[C:1](=[O:1])O.[N:3]>>[C:1](=[O:2])[N:3]");
		assertNotNull(rxn);
		validationResults = rxn.validateReaction();
		assertEquals(1, validationResults.getFirst());
		assertEquals(1, validationResults.getSecond());

		rxn = ChemicalReaction.ReactionFromSmarts("[C:1](=[O:2])[O:4].[N:3]>>[C:1](=[O:2])[N:3]");
		assertNotNull(rxn);
	}
	@Test (expected=ChemicalReactionException.class)
	public void test6Exceptions_1() {
		ChemicalReaction rxn = ChemicalReaction.ReactionFromSmarts("[C:1]Cl>>[C:1]");
		assertNotNull(rxn);
		rxn.runReactants(new ROMol_Vect());
	}
	@Test (expected=ChemicalReactionException.class)
	public void test6Exceptions_2() {
		ChemicalReaction rxn = ChemicalReaction.ReactionFromSmarts("[C:1]Cl>>[C:1]");
		assertNotNull(rxn);
		ROMol_Vect reactants = new ROMol_Vect();
		reactants.add(RWMol.MolFromSmiles("CC"));
		reactants.add(RWMol.MolFromSmiles("C"));
		rxn.runReactants(reactants);
	}
	@Test
	public void test6Exceptions_3() {
		ChemicalReaction rxn = ChemicalReaction.ReactionFromSmarts("[C:1]Cl>>[C:1]");
		assertNotNull(rxn);
		ROMol_Vect reactants = new ROMol_Vect();
		reactants.add(RWMol.MolFromSmiles("CCCl"));
		ROMol_Vect_Vect products = rxn.runReactants(reactants);
		assertEquals( 1,products.size() );
		assertEquals( 1,products.get(0).size() );
	}
	@Test
	public void test8Properties() {
		ChemicalReaction rxn = ChemicalReaction.ReactionFromSmarts("[O:1]>>[O:1][3#0]");
		assertNotNull(rxn );;
		ROMol_Vect reactants = new ROMol_Vect();
		reactants.add(RWMol.MolFromSmiles("CO"));
		ROMol_Vect_Vect ps = rxn.runReactants(reactants);
		assertEquals( 1,ps.size() );
		assertEquals( 1,ps.get(0).size() );
		// Can sanitize only RWMol objects
		RWMol p00 = new RWMol(ps.get(0).get(0));
		RDKFuncs.sanitizeMol(p00);
		assertEquals( 3,p00.getAtomWithIdx(1).getIsotope());
	}
	@Test
	public void test9AromaticityTransfer () {
		ROMol mol = RWMol.MolFromSmiles("c1ccc(C2C3(Cc4c(cccc4)C2)CCCC3)cc1");
		ChemicalReaction rxn = 
			ChemicalReaction.ReactionFromSmarts("[A:1]1~[*:2]~[*:3]~[*:4]~[*:5]~[A:6]-;@1>>[*:1]~[*:2]~[*:3]~[*:4]~[*:5]~[*:6]");
		ROMol_Vect reactants = new ROMol_Vect();
		reactants.add(mol);
		ROMol_Vect_Vect products = rxn.runReactants(reactants);
		assertEquals( 6,products.size() );

		for (int prodSet = 0; prodSet < products.size(); prodSet++) {
			assertEquals( 1,products.get(prodSet).size() );
			RWMol prod = new RWMol(products.get(prodSet).get(0));
			prod.sanitizeMol();
		}
	}

	@Test
	public void test10DotSeparation() {
		ROMol mol = RWMol.MolFromSmiles("C1ON1");
		ChemicalReaction rxn = 
			ChemicalReaction.ReactionFromSmarts("[C:1]1[O:2][N:3]1>>([C:1]1[O:2].[N:3]1)");
		ROMol_Vect reactants = new ROMol_Vect();
		reactants.add(mol);
		ROMol_Vect_Vect products = rxn.runReactants(reactants);
		assertEquals( 1,products.size() );
		for (int prodSet = 0; prodSet < products.size(); prodSet++) {
			assertEquals( 1,products.get(prodSet).size() );
			assertEquals( 3,products.get(prodSet).get(0).getNumAtoms() );
			assertEquals( 2,products.get(prodSet).get(0).getNumBonds() );
		}
	}

	@Test
	public void test11ImplicitProperties() {
		ROMol mol = RWMol.MolFromSmiles("CCO");
		ROMol mol2 = RWMol.MolFromSmiles("C[CH-]O");
		ChemicalReaction rxn = 
			ChemicalReaction.ReactionFromSmarts("[C:1]O>>[C:1]");
		ROMol_Vect reactants = new ROMol_Vect();
		reactants.add(mol);
		ROMol_Vect_Vect products = rxn.runReactants(reactants);
		assertEquals( 1,products.size() );
		for (int prodSet = 0; prodSet < products.size(); prodSet++) {
			assertEquals( 1,products.get(prodSet).size() );
			assertEquals("CC",products.get(prodSet).get(0).MolToSmiles());
		}
		reactants.clear();
		reactants.add(mol2);
		products = rxn.runReactants(reactants);
		assertEquals( 1,products.size() );
		for (int prodSet = 0; prodSet < products.size(); prodSet++) {
			assertEquals( 1,products.get(prodSet).size() );
			assertEquals("[CH2-]C",products.get(prodSet).get(0).MolToSmiles());
		}

		rxn.setImplicitPropertiesFlag(false);
		products = rxn.runReactants(reactants);
		assertEquals( 1,products.size() );
		for (int prodSet = 0; prodSet < products.size(); prodSet++) {
			assertEquals( 1,products.get(prodSet).size() );
			assertEquals("CC",products.get(prodSet).get(0).MolToSmiles());
		}
		reactants.clear();
		reactants.add(mol2);
		products = rxn.runReactants(reactants);
		assertEquals( 1,products.size() );
		for (int prodSet = 0; prodSet < products.size(); prodSet++) {
			assertEquals( 1,products.get(prodSet).size() );
			assertEquals("CC",products.get(prodSet).get(0).MolToSmiles());
		}
	}


    // @Test
	public void test99MemoryLeak() {
		ChemicalReaction rxn = 
			ChemicalReaction.ReactionFromSmarts("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]");
		assertNotNull(rxn );;
		assertEquals( 2,rxn.getNumReactantTemplates() );
		assertEquals( 1,rxn.getNumProductTemplates() );
		assertTrue(rxn.getImplicitPropertiesFlag());

		ROMol_Vect reacts = new ROMol_Vect();
		for (String react : new String[] {"C(=O)O","N"})
			reacts.add(RWMol.MolFromSmiles(react));		
		// Do not need the initReactantMatchers call here
                for(Integer i=0;i<1000000;i++){
                    ROMol_Vect_Vect prods = rxn.runReactants(reacts);;
                    assertEquals( 1,prods.size() );
                    assertEquals( 1,prods.get(0).size() );
                    assertEquals( 3,prods.get(0).get(0).getNumAtoms() );
                    //prods.get(0).get(0).delete();
                    //prods.get(0).delete();
                    //prods.delete();
                }

	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.ChemReactionTests");
	}

}
