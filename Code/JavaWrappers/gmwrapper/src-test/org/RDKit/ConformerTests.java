/* 
 * $Id: ConformerTests.java 131 2011-01-20 22:01:29Z ebakke $
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

public class ConformerTests extends GraphMolTest {

	@Test
	public void test0Conformers() {
		ROMol mol = RWMol.MolFromSmiles("CC");
		Conformer conf = new Conformer(2);
		conf.setAtomPos(0, new Point3D(-0.5, 0.0, 0.0));
		conf.setAtomPos(1, new Point3D(1.0, 0.0, 0.0));
		conf.setId(0);
		long cid = mol.addConformer(conf);

		assertEquals(0, cid);

		Conformer conf2 = mol.getConformer(0);
		assertEquals(cid, conf2.getId());
		Point3D pt1 = conf2.getAtomPos(0);
		assertPointEquals(pt1, new Point3D(-0.5, 0.0, 0.0));

		Point3D pt2 = conf2.getAtomPos(1);
		assertPointEquals(pt2, new Point3D(1.0, 0.0, 0.0));

		// conf and conf2 now reference the same object
		conf.setAtomPos(1, new Point3D(2.0, 0.0, 0.0));
		pt2 = conf2.getAtomPos(1);
		assertEquals(2.0, pt2.getX(), defaultDoubleTol);

		conf = new Conformer(2);
		conf.setAtomPos(0, new Point3D(-0.5, 0.0, 0.0));
		conf.setAtomPos(1, new Point3D(1.0, 0.0, 0.0));
		conf.setId(2);

		cid = mol.addConformer(conf, false);
		assertEquals(2, cid);
	}

	@Test
	public void test0AddHds() {
		ROMol mol = RWMol.MolFromSmiles("CC");
		Conformer conf = new Conformer(1);
		conf.setAtomPos(0, new Point3D(-0.5, 0.0, 0.0));
		conf.setAtomPos(1, new Point3D(1.0, 0.0, 0.0));
		long cid = mol.addConformer(conf);

		Conformer conf2 = mol.getConformer((int) cid);
		assertEquals(2, conf2.getNumAtoms());

		ROMol nmol = mol.addHs(false, true);
		Conformer conf3 = nmol.getConformer();
		assertEquals(8, conf3.getNumAtoms());
		assertEquals(2, conf2.getNumAtoms());

		double[][] targetCoords = { { -0.5, 0.0, 0.0 }, { 1.0, 0.0, 0.0 },
				{ -0.8667, 0.0, 1.03709 }, { -0.8667, 0.8981, -0.5185 },
				{ -0.8667, -0.8981, -0.5185 }, { 1.3667, 0.0, -1.0371 },
				{ 1.36667, 0.8981, 0.5185 }, { 1.36667, -0.8981, 0.5185 } };
		for (int i = 0; i < targetCoords.length; i++) {
			Point3D pt = conf3.getAtomPos(i);
			assertPointEquals(pt, new Point3D(targetCoords[i][0], targetCoords[i][1],
					targetCoords[i][2]));
		}
	}

	@Test
	public void test2Issue217() {
		String smi = "c1ccccc1";
		RWMol m = RWMol.MolFromSmiles(smi);
		m.addConformer(new Conformer(m.getNumAtoms()));
		assertEquals(1,m.getNumConformers());
		// We don't want an exception here
		m.MolToMolBlock();
	}
	@Test (expected=GenericRDKitException.class) 
	public void test3Exceptions() {
		String smi = "c1ccccc1";
		RWMol m = RWMol.MolFromSmiles(smi);
		m.addConformer(new Conformer());
		assertEquals(1,m.getNumConformers());
		// We just want the exception to be thrown
		m.getConformer(2);
	}

	@Test
	public void test4ConfTuple() {
		String smi = "c1ccccc1";
		RWMol m = RWMol.MolFromSmiles(smi);
		for (int i = 0; i < 10; i++) {
			Conformer c = new Conformer(m.getNumAtoms());
			c.setId(i);
			m.addConformer(c, false);
		}
		assertEquals(10, m.getNumConformers());

		for (int i = 0; i < m.getNumConformers(); i++)
			for (int j = 0; j < m.getConformer(i).getNumAtoms(); j++) {
				Point3D pt = m.getConformer(i).getAtomPos(j);
				assertPointEquals(pt, new Point3D(0.0, 0.0, 0.0));
			}
		for (int i = 0; i < 10; i++)
			m.removeConformer(i);
		assertEquals(0, m.getNumConformers());
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.ConformerTests");
	}
}
