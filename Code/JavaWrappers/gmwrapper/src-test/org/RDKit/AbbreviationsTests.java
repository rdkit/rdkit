/* 
 *
 *  Copyright (c) 2019 Greg Landrum and T5 Informatics GmbH
 *  All rights reserved.
 * 
 *  This file is part of the RDKit.
 *  The contents are covered by the terms of the BSD license
 *  which is included in the file license.txt, found at the root
 *  of the RDKit source tree.
 */
package org.RDKit;

import static org.junit.Assert.*;
import org.junit.*;
    
public class AbbreviationsTests extends GraphMolTest {

	@Before public void setUp() {
	}

	@Test 
	public void test1Basics() {
        AbbreviationDefinition_Vect abbrevs = RDKFuncs.getDefaultAbbreviations();
        RWMol mol = RWMol.MolFromSmiles("C1CCC1C(F)(F)F");
        assertEquals(mol.getNumAtoms(),8);

		RDKFuncs.condenseMolAbbreviations(mol,abbrevs);
        // no changes here due to the threshold
        assertEquals(mol.getNumAtoms(),8);

		RDKFuncs.condenseMolAbbreviations(mol,abbrevs, 1.0);
        assertEquals(mol.getNumAtoms(),5);
        assertEquals(RDKFuncs.MolToCXSmiles(mol),"*C1CCC1 |$CF3;;;;$|");	  
	}
  
	@Test 
	public void test2LinkerBasics() {
        AbbreviationDefinition_Vect abbrevs = RDKFuncs.getDefaultLinkers();
        RWMol mol = RWMol.MolFromSmiles("COCCOCCOCCOCCCl");
        assertEquals(mol.getNumAtoms(),14);

		RDKFuncs.condenseMolAbbreviations(mol,abbrevs);
        // no changes here due to the threshold
        assertEquals(mol.getNumAtoms(),14);

		RDKFuncs.condenseMolAbbreviations(mol,abbrevs, 1.0);
        assertEquals(mol.getNumAtoms(),3);
        assertEquals(RDKFuncs.MolToCXSmiles(mol),"C*Cl |$;PEG4;$|");
		  
	}

	@Test 
	public void test3Matching() {
        AbbreviationDefinition_Vect abbrevs = RDKFuncs.getDefaultAbbreviations();
        RWMol mol = RWMol.MolFromSmiles("C1CCC1C(F)(F)F");
        assertEquals(mol.getNumAtoms(),8);

		AbbreviationMatch_Vect matches = RDKFuncs.findApplicableAbbreviationMatches(mol,abbrevs,1.0);
        assertEquals(matches.size(),1);
        assertEquals(matches.get(0).getAbbrev().getLabel(),"CF3");

        RDKFuncs.applyMatches(mol,matches);
        assertEquals(mol.getNumAtoms(),5);
        assertEquals(RDKFuncs.MolToCXSmiles(mol),"*C1CCC1 |$CF3;;;;$|");	  

	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.AbbreviationsTests");
	}

}
