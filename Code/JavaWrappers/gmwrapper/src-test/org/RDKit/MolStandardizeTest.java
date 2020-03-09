package org.RDKit;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.*;

/**
 * JUnit tests for MolStandardize wrappers.
 */
 public class MolStandardizeTest extends GraphMolTest {
	@Test
	public void testStandardize1() {
		assertEquals("fail", RDKFuncs.standardizeSmiles("[Na]OC(=O)c1ccccc1"),"O=C([O-])c1ccccc1.[Na+]");
	}
	public static void main(String args[]) {
						org.junit.runner.JUnitCore.main("org.RDKit.MolStandardizeTest");
	}
}
