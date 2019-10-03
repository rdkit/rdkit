package org.RDKit;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.*;

/**
 * JUnit tests for MolHash wrappers.
 */
 public class MolHashTest extends GraphMolTest {
	@Test
	public void testMolHash1() {
		ROMol m1;
        m1 = RWMol.MolFromSmiles("C1CCCC(O)C1c1ccnc(OC)c1");
		// we just test a few to make sure that the enum works
		assertEquals("***1****(*2*****2*)*1", RDKFuncs.MolHash(new RWMol(m1),HashFunction.AnonymousGraph));
		assertEquals("COc1cc(C2CCCCC2O)ccn1", RDKFuncs.MolHash(new RWMol(m1),HashFunction.CanonicalSmiles));
		assertEquals("COC1CC(C2CCCCC2O)CCN1", RDKFuncs.MolHash(new RWMol(m1),HashFunction.ElementGraph));
	}
	public static void main(String args[]) {
						org.junit.runner.JUnitCore.main("org.RDKit.MolHashTest");
	}
}
