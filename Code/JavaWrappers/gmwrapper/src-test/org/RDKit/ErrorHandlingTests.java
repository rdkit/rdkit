/* 
 * $Id: ErrorHandlingTests.java 131 2011-01-20 22:01:29Z ebakke $
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

// Uncomment the following line and enable the creation of the
// ErrorGenerator class in  .../gmwrapper/CMakeLists.txt to
// test behavior around low-level C++ errors.  Also uncomment the
// appropriate tests below.
//import org.RDKit.ErrorGenerator;

public class ErrorHandlingTests extends GraphMolTest {

	@Test
	public void testDataGetSetFailure() {
		ROMol m = RWMol.MolFromSmiles("CCOC");
		// m.setProp("monkey", "Spider");
		try {
			m.getProp("monkey");
			fail("Should not be able to access 'monkey'");
		}
		catch (Exception e)
		{
			assertEquals("class org.RDKit.KeyErrorException",
					e.getClass().toString());
			String what = ((KeyErrorException) e).key();
			assertEquals("monkey", what);
			what = ((KeyErrorException) e).message();
			assertEquals("Unknown key: monkey", what);
		}
	}

	@Test(expected=GenericRDKitException.class)
	public void testAtomIndexInMolecule () {
		ROMol mol = RWMol.MolFromSmiles("CCO");
		Atom a = mol.getAtomWithIdx(mol.getNumAtoms());
		@SuppressWarnings("unused")
		int n = a.getAtomicNum();
	}

	@Test
	public void testAtomIndexInMoleculeCaught () {
		ROMol mol = RWMol.MolFromSmiles("CCO");
		try {
			@SuppressWarnings("unused")
			Atom a = mol.getAtomWithIdx(mol.getNumAtoms());
			fail();
		}
		catch (GenericRDKitException e)
		{
			String what = e.message();
			assertEquals("Unknown exception", what);
			String where = e.getStackTrace()[0].getMethodName();
			assertEquals("ROMol_getAtomWithIdx", where);
		}
	}

	// The following tests require the ErrorGenerator class
	// described above.  Uncomment them when that class is available.
	// Note that several tests are marked @Ignore because they
	// demonstrate situations which can't be trapped and cause
	// the JVM to crash.
	/*
	@Test(expected=GenericRDKitException.class)
	public void testBadAlloc_1() {
		ErrorGenerator eg = new ErrorGenerator();
		eg.badAlloc_1();
	}

	@Test
	public void testBadAlloc_1b() {
		ErrorGenerator eg = new ErrorGenerator();
		try {
			eg.badAlloc_1();
			fail("Should not be able to allocate");
		}
		catch (GenericRDKitException e)
		{
			assertEquals("class org.RDKit.GenericRDKitException",
					e.getClass().toString());
			String msg = e.message();
			assertEquals("Unknown exception", msg);
			StackTraceElement[] st = e.getStackTrace();
			assertEquals("ErrorGenerator_badAlloc_1", 
					st[0].getMethodName());
		}

	}

	@Test(expected=GenericRDKitException.class)
	public void testBadAlloc_2() {
		ErrorGenerator eg = new ErrorGenerator();
		eg.badAlloc_2();
	}
	// Here is an example of errors that can't be caught

	@Test(expected=GenericRDKitException.class)
	@Ignore // Note that this *WILL* cause a seg fault and a JVM crash
	// -- can't help it
	public void testBadCall() {
		ErrorGenerator eg = new ErrorGenerator();
		eg.badCall();
	}

	@Test(expected=GenericRDKitException.class)
	@Ignore // Note that this *WILL* cause a seg fault and a JVM crash
	// -- can't help it
	public void testBadAccess() {
		ErrorGenerator eg = new ErrorGenerator();
		int i = eg.badAccess();
		// Force a reference to i so that the compiler doesn't
		// optimize away the above call.
		@SuppressWarnings("unused")
		int j = i + 1;

	}
	*/

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.ErrorHandlingTests");
	}
}
