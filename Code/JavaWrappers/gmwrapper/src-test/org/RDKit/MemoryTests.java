/* 
 * $Id: MemoryTests.java 131 2011-01-20 22:01:29Z ebakke $
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

/**
 * A few tests to verify memory-relates issues between Java and C++.
 * 
 * Verify garbage collection correctly frees Java and C++ memory.
 * Call it manually, otherwise it only gets called at the end of all tests
 */
public class MemoryTests extends GraphMolTest {

	Runtime r = Runtime.getRuntime();
	static final int ITERATIONS = 50000; 

	@SuppressWarnings("unused")
	@Test
	public void testGarbageCollection1() {
		for(int i=0;i<ITERATIONS;i++){
			ROMol m1 = RWMol.MolFromSmiles("C1=CC=CC=C1");
			if((i%5000)==0){
				r.gc();
				r.runFinalization();
			}
		}
	}
	
	@SuppressWarnings("unused")
	@Test 
	public void testGarbageCollection2(){
		int i = 0;
		ROMol mol;
		String smi = "C1CC1";
		for(i=0;i<ITERATIONS;i++){
			mol = RWMol.MolFromSmiles(smi,0,true);
			if((i%5000)==0){
				r.gc();
				r.runFinalization();
			}
		}
	}

	@Test
	public void testConformerOwnership() {
		for (int i = 0; i < 50000; i++) {
			ROMol mol = RWMol.MolFromSmiles("CC");
			
			for(int j=0; j<100;j++) {
			Conformer conf = new Conformer(2);
			conf.setAtomPos(0, new Point3D(-0.5, 0.0, 0.0));
			conf.setAtomPos(1, new Point3D(1.0, 0.0, 0.0));
			conf.setId(j);
			mol.addConformer(conf);

			conf = new Conformer(2);
			conf.setAtomPos(0, new Point3D(-0.5, 0.0, 0.0));
			conf.setAtomPos(1, new Point3D(1.0, 0.0, 0.0));
			conf.setId(j+1);
			mol.addConformer(conf, false);
			}
			if (i % 1000 == 0) {
				r.gc();
				r.runFinalization();
			}
		}
	}

	@Test 
	public void testAtomsAfterGarbageCollection() {
		String smiles="c1ccccc1";
		RWMol mol1 = RWMol.MolFromSmiles(smiles);

		assertEquals(6, mol1.getAtomWithIdx(0).getAtomicNum());
		assertNotNull(mol1.getAtomWithIdx(1));
		for( int i = 0; i< mol1.getNumAtoms(); i++) {
			assertNotNull( mol1.getAtomWithIdx(i));
			assertFalse( 47 == mol1.getAtomWithIdx(i).getAtomicNum() );

			Atom atom = new Atom(47);
			mol1.replaceAtom(i,atom,false);
			assertEquals(47, mol1.getAtomWithIdx(i).getAtomicNum());
		}
		
		// Atom's now out of scope and eligible for garbage collection!
		// Make sure we still have access to it, even after the Java object was destroyed.
		Runtime.getRuntime().gc();
		Runtime.getRuntime().runFinalization();
		for( int i = 0; i < mol1.getNumAtoms(); i++) {
			Atom atom = mol1.getAtomWithIdx(i);
			assertNotNull(atom);
			assertEquals( 47 , atom.getAtomicNum() );
		}
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.MemoryTests");
	}
}


