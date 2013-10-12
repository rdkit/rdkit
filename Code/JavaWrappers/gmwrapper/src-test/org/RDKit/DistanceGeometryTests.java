/* 
 * $Id: DistanceGeometryTests.java 131 2011-01-20 22:01:29Z ebakke $
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

import java.io.File;

import org.junit.*;

public class DistanceGeometryTests extends GraphMolTest {

	@Before
	public void setUp() {
	}
	@Test
	public void testBasicEmbedding() {
		String[] smiString = new String[] {
				"CC1=C(C(C)=CC=C2)C2=CC=C1", "c1ccccc1C", "C/C=C/CC",
				"C/12=C(\\CSC2)Nc3cc(n[n]3C1=O)c4ccccc4", "C1CCCCS1(=O)(=O)", "c1ccccc1",
				"C1CCCC1", "C1CCCCC1",
				"C1CC1(C)C", "C12(C)CC1CC2"};
		String fname = new File(getRdBase(), 
		"Code/GraphMol/DistGeomHelpers/test_data/initCoords.sdf").getPath();
		SDMolSupplier sdsup = new SDMolSupplier(fname);

		for (String smi : smiString) {
			RWMol m = RWMol.MolFromSmiles(smi, 0, true); 
			DistanceGeom.EmbedMolecule(m, 10, 1,true,false,2,true,1,null,1e-2);

			ROMol m2 = sdsup.next();
			if(m2 != null){
				int nat = (int) m.getNumAtoms();

				Conformer conf1 = m.getConformer(0);
				Conformer conf2 = m2.getConformer(0);
				for (int i = 0; i < nat; i++) {
					Point3D pt1i = conf1.getAtomPos(i);
					Point3D pt2i = conf2.getAtomPos(i);
					for(int j=i+1;j<nat;j++){
						Point3D pt1j = conf1.getAtomPos(j);
						Point3D pt2j = conf2.getAtomPos(j);
						double d1 = pt1j.minus(pt1i).length();
						double d2 = pt2j.minus(pt2i).length();
						if(m.getBondBetweenAtoms(i,j) != null){
							assertEquals(d1, d2, 0.06);
                                                } else {
							assertEquals(d1, d2, 0.15);
                                                }
					}
				}
			}
		}
	}

	private void computeDistMat(Point3D_Vect origCoords, DoubleSymmMatrix distMat) {
		int N = (int) origCoords.size();
		assertEquals(N, distMat.numRows());
		int i, j;
		Point3D pti, ptj;
		double d;
		for (i = 1; i < N; i++) {
			pti = origCoords.get(i);
			for (j = 0; j < i; j++) {
				ptj = origCoords.get(j);
				ptj = ptj.minus(pti);
				d = ptj.length();
				distMat.setVal(i,j, d);
			}
		}
	}

	private void computeMolDmat(ROMol mol, DoubleSymmMatrix distMat) {
		Point3D_Vect origCoords = new Point3D_Vect();
		int nat = (int) mol.getNumAtoms();
		Conformer conf = mol.getConformer(0);
		for (int i = 0; i < nat; i++) {
			origCoords.add(conf.getAtomPos(i));
		}
		computeDistMat(origCoords, distMat);
	}

	@Test
	public void testInRingDistances() {
		String smi = "Cc1c(C=CC(C)N2)c2[nH]n1";
		ROMol mol = RWMol.MolFromSmiles(smi, 0, true);
		int cid;
		long nat = mol.getNumAtoms();
		BoundsMatrix bm = new BoundsMatrix(nat);
		RDKFuncs.initBoundsMat(bm, 0.0, 1000.0);
		DistanceGeom.SetTopolBounds(mol, bm);
		cid = DistanceGeom.EmbedMolecule(mol, 10, 1);
		assertTrue(cid>-1);
		DoubleSymmMatrix dmat = new DoubleSymmMatrix(nat, 0.0);
		computeMolDmat(mol, dmat);

		assertTrue((bm.getUpperBound(0,9) - bm.getLowerBound(0,9)) < 0.13);
		assertTrue((bm.getUpperBound(0,9) - dmat.getVal(0,9) > -0.1 )
				&& (bm.getLowerBound(0,9) - dmat.getVal(0,9) < 0.10));

		assertTrue((bm.getUpperBound(10,7) - bm.getLowerBound(10,7)) < 0.13);
		assertTrue((bm.getUpperBound(10,7) - dmat.getVal(10,7) > -0.1 )
				&& (bm.getLowerBound(10,7) - dmat.getVal(10,7) < 0.10 ));

		assertTrue((bm.getUpperBound(2,5) - bm.getLowerBound(2,5)) < 0.20);
		assertTrue((bm.getUpperBound(2,5) - dmat.getVal(2,5) > -0.1 )
				&& (bm.getLowerBound(2,5) - dmat.getVal(2,5) < 0.10 ));

		assertTrue((bm.getUpperBound(8,4) - bm.getLowerBound(8,4)) > 1.);
		assertTrue((bm.getUpperBound(8,4) - bm.getLowerBound(8,4)) < 1.2);
		assertTrue((bm.getUpperBound(8,4) - dmat.getVal(8,4) > -0.1 )
				&& (bm.getLowerBound(8,4) - dmat.getVal(8,4) < 0.10));

		assertTrue((bm.getUpperBound(8,6) - bm.getLowerBound(8,6)) > 1.0);
		assertTrue((bm.getUpperBound(8,6) - bm.getLowerBound(8,6)) < 1.2);
		assertTrue((bm.getUpperBound(8,6) - dmat.getVal(8,6) > -0.1)
				&& (bm.getLowerBound(8,6) - dmat.getVal(8,6) < 0.10 ));


	}
	@Test
	public void testChainOfSingles_1() {
		// chain of singles
		String smi = "CCCC";
		ROMol mol = RWMol.MolFromSmiles(smi, 0, true);
		long nat = mol.getNumAtoms();
		BoundsMatrix bm = new BoundsMatrix(nat);
		RDKFuncs.initBoundsMat(bm, 0.0, 1000.0);
		DistanceGeom.SetTopolBounds(mol, bm);
		int cid = DistanceGeom.EmbedMolecule(mol, 10, 1);
		assertTrue(cid>-1);
		DoubleSymmMatrix dmat = new DoubleSymmMatrix(nat, 0.0);
		computeMolDmat(mol, dmat);
		assertTrue( (bm.getUpperBound(0,3) - bm.getLowerBound(0,3)) > 1.0);
		assertTrue( (bm.getUpperBound(0,3) - bm.getLowerBound(0,3)) < 1.3);
		assertTrue( (bm.getUpperBound(0,3) - dmat.getVal(0,3) > -0.1)
				&& (bm.getLowerBound(0,3) - dmat.getVal(0,3) < 0.10 ));
	}
	@Test
	public void testCompareEmbeddingToCombiCoords () {
		String fname = new File(getRdBase(), 
		"Code/GraphMol/DistGeomHelpers/test_data/combi_coords.sdf").getPath();
		SDMolSupplier sdsup = new SDMolSupplier(fname);

		while (!sdsup.atEnd()) {
			ROMol mol = sdsup.next();
			Point3D_Vect origCoords = new Point3D_Vect(); 
			long nat = mol.getNumAtoms();
			Conformer conf = mol.getConformer(0);
			for (int i = 0; i < nat; i++) {
				origCoords.add(conf.getAtomPos(i));
			}
			DoubleSymmMatrix distMat = new DoubleSymmMatrix(nat, 0.0);
			computeDistMat(origCoords, distMat);
			assertTrue(DistanceGeom.ComputeInitialCoords(distMat, origCoords));

			DoubleSymmMatrix distMatNew = new DoubleSymmMatrix(nat, 0.0);
			computeDistMat(origCoords, distMatNew);

			for (int i = 1; i < nat; i++) {
				for (int j = 0; j < i; j++) {
					assertEquals(distMat.getVal(i,j), distMatNew.getVal(i,j), 0.01);
				}
			}
		}			
	}
	@Test
	public void testIssue215 () {
		ROMol m = RWMol.MolFromSmiles("C=C1C2CC1C2");
		assertNotNull(m);
		BoundsMatrix bm = new BoundsMatrix(m.getNumAtoms());
		assertNotNull(bm);
		RDKFuncs.initBoundsMat(bm,0.0,1000.0);
		DistanceGeom.SetTopolBounds(m, bm);

		// this was the specific problem:
		assertTrue(bm.getUpperBound(0,4)<100.0);

		assertTrue(RDKFuncs.triangleSmoothBounds(bm));

	}

	@Test
	public void test15Dists() {

		ROMol m = RWMol.MolFromSmiles("c1ccccc1C");
		int nat = (int) m.getNumAtoms();
		BoundsMatrix mat = new BoundsMatrix(nat);
		RDKFuncs.initBoundsMat(mat);
		DistanceGeom.SetTopolBounds(m, mat);
		assertEquals(mat.getUpperBound(2,6), 4.32, 0.01);
		assertEquals(mat.getLowerBound(2,6), 4.16, 0.01);

		m = RWMol.MolFromSmiles("CC1=C(C(C)=CC=C2)C2=CC=C1");
		nat = (int) m.getNumAtoms();
		mat = new BoundsMatrix(nat);
		RDKFuncs.initBoundsMat(mat);
		DistanceGeom.SetTopolBounds(m, mat);

		assertEquals(mat.getLowerBound(0,4), 2.31, 0.01);
		assertEquals(mat.getUpperBound(0,4), 2.47, 0.01);
		assertEquals(mat.getLowerBound(4,11), 4.11, 0.01);
		assertEquals(mat.getUpperBound(4,11), 4.27, 0.01) ;


		m = RWMol.MolFromSmiles("C/C=C/C=C/C", 0, true);
		nat = (int) m.getNumAtoms();

		mat = new BoundsMatrix(nat);
		RDKFuncs.initBoundsMat(mat);
		DistanceGeom.SetTopolBounds(m, mat);

		assertEquals(mat.getLowerBound(0,4), 4.1874, 0.01);
		assertEquals(mat.getUpperBound(0,4), 4.924, 0.01);
		assertEquals(mat.getLowerBound(1,5), 4.1874, 0.01);
		assertEquals(mat.getUpperBound(1,5), 4.924, 0.01) ;
	}
	@Test
	public void testMultipleConfs() {
		String smi = "CC(C)(C)c(cc1)ccc1c(cc23)n[n]3C(=O)/C(=C\\N2)C(=O)OCC";
		ROMol m = RWMol.MolFromSmiles(smi, 0, true);
		Int_Vect cids = DistanceGeom.EmbedMultipleConfs(m, 10, 30, 100, true, false,-1);
		SDWriter writer = new SDWriter("junk.sdf");
		double energy;

		for (int i = 0; i < cids.size(); i++ ) {
			writer.write(m, cids.get(i));
			ForceField ff = ForceField.UFFGetMoleculeForceField(m, 10, cids.get(i));
			ff.initialize();
			energy = ff.calcEnergy();
			assertTrue(energy>100.0);
			assertTrue(energy<300.0);
		}

	}
	@Test
	public void testConstrainedEmbedding() {
		String fname = new File(getRdBase(), 
		"Code/GraphMol/DistGeomHelpers/test_data/constrain1.sdf").getPath();
		SDMolSupplier sdsup = new SDMolSupplier(fname);

		ROMol ref=sdsup.next();
		ROMol test = new ROMol(ref);
		Int_Point3D_Map coords = new Int_Point3D_Map();
		coords.set(0, ref.getConformer().getAtomPos(0));
		coords.set(1, ref.getConformer().getAtomPos(1));
		coords.set(2, ref.getConformer().getAtomPos(2));
		coords.set(3, ref.getConformer().getAtomPos(3));
		coords.set(4, ref.getConformer().getAtomPos(4));

		int cid = DistanceGeom.EmbedMolecule(test, 30,23,true,false,2.,true,1, coords);
		assertTrue(cid>-1);

		Match_Vect alignMap = new Match_Vect();
		alignMap.add(new Int_Pair(0, 0));
		alignMap.add(new Int_Pair(1, 1));
		alignMap.add(new Int_Pair(2, 2));
		alignMap.add(new Int_Pair(3, 3));
		alignMap.add(new Int_Pair(4, 4));

		double ssd = test.alignMol(ref,-1,-1,alignMap);
		assertTrue(ssd<0.1);

		test = sdsup.next();
		coords.clear();
		coords.set(4, ref.getConformer().getAtomPos(0));
		coords.set(5, ref.getConformer().getAtomPos(1));
		coords.set(6, ref.getConformer().getAtomPos(2));
		coords.set(7, ref.getConformer().getAtomPos(3));
		coords.set(8, ref.getConformer().getAtomPos(4));
		cid = DistanceGeom.EmbedMolecule(test, 30,23,true,false,2.,true,1,coords);
		assertTrue(cid>-1);

		alignMap.clear();
		alignMap.add(new Int_Pair(4,0));
		alignMap.add(new Int_Pair(5,1));
		alignMap.add(new Int_Pair(6,2));
		alignMap.add(new Int_Pair(7,3));
		alignMap.add(new Int_Pair(8,4));
		ssd = test.alignMol(ref,-1,-1,alignMap);
		assertTrue(ssd < 0.1);
	}


	@Test
	public void testIssue2835784() {
		String smi="C1C=C1";
		RWMol m = RWMol.MolFromSmiles(smi);
		int cid = DistanceGeom.EmbedMolecule(m);
		assertTrue(cid >= 0);
		Int_Vect cids=DistanceGeom.EmbedMultipleConfs(m, 10);
		assertEquals(cids.size(), 10);
		// Make sure there are no id values of -1
		for (int i = 0; i < cids.size(); i++)
			assertFalse("No id of -1 expected", -1 == cids.get(i));
	}

	// from testDistGeom.py
	@Test
	public void test7ConstrainedEmbedding() {
		File ofile = new File(getRdBase(), 
				"Code/GraphMol/DistGeomHelpers/test_data/constrain1.sdf");
		SDMolSupplier suppl = new SDMolSupplier(ofile.getPath());
        ROMol ref = suppl.next();
        ROMol probe = new ROMol(ref);

        Int_Point3D_Map cmap = new Int_Point3D_Map();
        for (int i = 0; i < 5; i++)
        	cmap.set(i, ref.getConformer().getAtomPos(i));
        int ci = DistanceGeom.EmbedMolecule(probe, 
        		0, 23, true, false, 2.0, false, 1, cmap);
        assertTrue(ci>-1);
        Match_Vect algMap = new Match_Vect();
        for (int i = 0; i < 5; i++)
        algMap.add(new Int_Pair(i, i));
        double ssd = probe.alignMol(ref, -1, -1, algMap);
        assertTrue(ssd<0.1);
	}
	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.DistanceGeometryTests");
	}

}
