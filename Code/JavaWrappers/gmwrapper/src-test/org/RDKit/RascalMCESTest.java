/* 
 * 
 *
 *  Copyright (c) 2024, Greg Landrum
 *  All rights reserved.
 * 
 *   @@ All Rights Reserved @@
 *  This file is part of the RDKit.
 *  The contents are covered by the terms of the BSD license
 *  which is included in the file license.txt, found at the root
 *  of the RDKit source tree.
 */
package org.RDKit;

import static org.junit.Assert.*;

import java.io.*;
import java.util.ArrayList;

import org.junit.*;

public class RascalMCESTest extends GraphMolTest {

	private File baseTestPath;

	@Before 
	public void setUp() {
		File base = getRdBase();
		baseTestPath = new File(base, "Contrib" + File.separator + "Fastcluster"+ File.separator + "cdk2.smi");
	}

	@After
	public void tearDown() {
	}

	@Test
	public void test1Rascal() {
		ROMol m1 = RWMol.MolFromSmiles("CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C");
		ROMol m2 = RWMol.MolFromSmiles("CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O");
        RascalOptions options = new RascalOptions();
        options.setSimilarityThreshold(0.6);
        RascalResult_Vect res = RDKFuncs.rascalMCES(m1, m2, options);
        assertEquals(res.size(),1);
        assertEquals(res.get(0).getSmarts(),"CC12CCC(-C(-C1CCC2O)-CC-[#6])-[#6]");
	}

	@Test
	public void test2RascalButina() {
		SmilesMolSupplier suppl = new SmilesMolSupplier(baseTestPath.getPath(),"\t", 1, 0, false);
        ROMol_Vect ms = new ROMol_Vect();
		do { 
			ms.add(suppl.next());
		}
		while (!suppl.atEnd());

        Unsigned_Vect_Vect res = RascalApp.RascalButinaCluster(ms);
        assertEquals(res.size(),29);
        assertEquals(res.get(0).size(),6);
        assertEquals(res.get(1).size(),6);

    }

	@Test
	public void test3RascalCluster() {
		SmilesMolSupplier suppl = new SmilesMolSupplier(baseTestPath.getPath(),"\t", 1, 0, false);
        ROMol_Vect ms = new ROMol_Vect();
		do { 
			ms.add(suppl.next());
		}
		while (!suppl.atEnd());

        Unsigned_Vect_Vect res = RascalApp.RascalCluster(ms);
        assertEquals(res.size(),8);
        assertEquals(res.get(0).size(),7);
        assertEquals(res.get(1).size(),7);

    }

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.RascalMCESTest");
	}

}
