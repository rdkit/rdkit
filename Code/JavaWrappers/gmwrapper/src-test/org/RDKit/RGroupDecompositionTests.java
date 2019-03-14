/* 
 *
 *  Copyright (c) 2019 Greg Landrum
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
    
public class RGroupDecompositionTests extends GraphMolTest {
	private ROMol mol;
    private ROMol m;

	@Before public void setUp() {
	}

	@Test 
	public void test1Basics() {
            ROMol core = RWMol.MolFromSmiles("c1ccccc1");
            RGroupDecomposition decomp = new RGroupDecomposition(core);

            m = RWMol.MolFromSmiles("c1(Cl)ccccc1");
            assertEquals(0,decomp.add(m));
            m = RWMol.MolFromSmiles("c1c(Cl)cccc1");
            assertEquals(1,decomp.add(m));
            assertTrue(decomp.process());

	}
  
	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.RGroupDecompositionTests");
	}

}
