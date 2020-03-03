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
    
public class ScaffoldNetworkTests extends GraphMolTest {

	@Before public void setUp() {
	}

	@Test 
	public void test1Basics() {
		ROMol_Vect ms = new ROMol_Vect();
		ms.add(RWMol.MolFromSmiles("c1ccccc1CC1NC(=O)CCC1"));
		ScaffoldNetworkParams ps = new ScaffoldNetworkParams();

		ScaffoldNetwork net = RDKFuncs.createScaffoldNetwork(ms,ps);

		assertEquals(net.getNodes().size(),9);
		assertEquals(net.getCounts().size(),net.getNodes().size());
		assertEquals(net.getEdges().size(),8);   
	}
  
	@Test 
	public void test2Basics() {
		ROMol_Vect ms = new ROMol_Vect();
		ms.add(RWMol.MolFromSmiles("c1ccccc1CC1NC(=O)CCC1"));
		ScaffoldNetworkParams ps = new ScaffoldNetworkParams();
		ps.setIncludeScaffoldsWithAttachments(false);

		ScaffoldNetwork net = RDKFuncs.createScaffoldNetwork(ms,ps);

		assertEquals(net.getNodes().size(),5);
		assertEquals(net.getCounts().size(),net.getNodes().size());
		assertEquals(net.getEdges().size(),4);   
	}
  
	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.ScaffoldNetworkTests");
	}

}
