/* 
 *  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
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

import org.junit.Test;

public class AlignTests extends GraphMolTest {


	@Test
	public void testO3ABasic () {
            String fname = new File(getRdBase(),
                                    "Code/GraphMol/MolAlign/test_data/ref_e2.sdf").getPath();
            SDMolSupplier sdsup = new SDMolSupplier(fname);
            ROMol m1 = sdsup.next();
            ROMol m2 = sdsup.next();
            Double_Pair res = m1.O3AAlignMol(m2);
            assertEquals(res.getFirst(),0.049,.001);
            assertEquals(res.getSecond(),119.98,.01);
            
	}
    
	@Test
	public void testgetAlignmentTransform(){
		//Test alignment with Transform base on GraphMol/MolAlign/testMolAlign.cpp#test1MolAlign()
		String fname0 = new File(getRdBase(),
					 "Code/GraphMol/MolAlign/test_data/1oir.mol").getPath();
		String fname1 = new File(getRdBase(),
					 "Code/GraphMol/MolAlign/test_data/1oir_conf.mol").getPath();
		ROMol m0 = RWMol.MolFromMolFile(fname0);
		ROMol m1 = RWMol.MolFromMolFile(fname1);
		Transform3D trans = new Transform3D();
		double res = m0.getAlignmentTransform(m1, trans);
		assertEquals(res, 1.0345, 0.001);
		m0.delete();
		m1.delete();
	}

	private String alignDataFile(String name) {
		return new File(getRdBase(),
				"Code/GraphMol/MolAlign/test_data/" + name).getPath();
	}

	private Match_Vect oirAtomMap() {
		int[][] pairs = {{18, 27}, {13, 23}, {21, 14}, {24, 7}, {9, 19}, {16, 30}};
		Match_Vect atomMap = new Match_Vect();
		for (int[] p : pairs) {
			atomMap.add(new Int_Pair(p[0], p[1]));
		}
		return atomMap;
	}

	// the tests below use the free functions from MolAlign/AlignMolecules.h
	// and MolAlign/O3AAlignMolecules.h and are based on
	// Code/GraphMol/MolAlign/Wrap/testMolAlign.py
	@Test
	public void testAlignMolFunction() {
		ROMol mol1 = RWMol.MolFromMolFile(alignDataFile("1oir.mol"));
		ROMol mol2 = RWMol.MolFromMolFile(alignDataFile("1oir_conf.mol"));
		double rmsd = RDKFuncs.alignMol(mol2, mol1);
		assertEquals(0.6578, rmsd, 0.0001);

		// the probe should now have the coordinates of the pre-aligned file
		ROMol mol3 = RWMol.MolFromMolFile(alignDataFile("1oir_trans.mol"));
		Conformer conf2 = mol2.getConformer();
		Conformer conf3 = mol3.getConformer();
		for (int i = 0; i < mol2.getNumAtoms(); ++i) {
			Point3D p2 = conf2.getAtomPos(i);
			Point3D p3 = conf3.getAtomPos(i);
			assertEquals(p3.getX(), p2.getX(), 0.001);
			assertEquals(p3.getY(), p2.getY(), 0.001);
			assertEquals(p3.getZ(), p2.getZ(), 0.001);
		}

		Transform3D trans = new Transform3D();
		rmsd = RDKFuncs.getAlignmentTransform(mol2, mol1, trans);
		assertEquals(0.6579, rmsd, 0.0001);
	}

	@Test
	public void testAlignMolAtomMapAndWeights() {
		ROMol mol1 = RWMol.MolFromMolFile(alignDataFile("1oir.mol"));
		ROMol mol2 = RWMol.MolFromMolFile(alignDataFile("1oir_conf.mol"));
		Match_Vect atomMap = oirAtomMap();
		double rmsd = RDKFuncs.alignMol(mol2, mol1, 0, 0, atomMap);
		assertEquals(0.8525, rmsd, 0.0001);

		mol2 = RWMol.MolFromMolFile(alignDataFile("1oir_conf.mol"));
		DoubleVector wts = new DoubleVector(6, 1.0);
		wts.setVal(5, 2.0);
		rmsd = RDKFuncs.alignMol(mol2, mol1, 0, 0, atomMap, wts);
		assertEquals(0.9513, rmsd, 0.0001);
	}

	@Test
	public void testAlignMolConformers() {
		ROMol mol = RWMol.MolFromSmiles("C1CC1CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4");
		Int_Vect cids = DistanceGeom.EmbedMultipleConfs(mol, 10, 30, 100);
		assertEquals(10, cids.size());
		UInt_Vect aids = new UInt_Vect();
		for (int aid : new int[] {12, 13, 14, 15, 16, 17, 18}) {
			aids.add(aid);
		}
		Double_Vect rmsVals = new Double_Vect();
		RDKFuncs.alignMolConformers(mol, aids, null, null, false, 50, rmsVals);
		assertEquals(mol.getNumConformers() - 1, rmsVals.size());
		for (int i = 0; i < rmsVals.size(); ++i) {
			assertTrue(rmsVals.get(i) >= 0.0);
		}
	}

	@Test
	public void testBestRMS() {
		SDMolSupplier suppl = new SDMolSupplier(alignDataFile("probe_mol.sdf"), true, false);
		suppl.moveTo(1);
		ROMol prb = suppl.next();
		ROMol ref = suppl.next();
		ROMol prbCopy = new ROMol(prb);

		double rmsdInPlace = RDKFuncs.CalcRMS(new ROMol(prb), ref);
		assertEquals(2.6026, rmsdInPlace, 0.001);
		// alignMol() would give 2.50561 here
		double rmsd = RDKFuncs.getBestRMS(prb, ref);
		assertEquals(2.43449, rmsd, 0.001);

		Transform3D bestTrans = new Transform3D();
		Match_Vect bestMatch = new Match_Vect();
		double rmsdCopy = RDKFuncs.getBestAlignmentTransform(prbCopy, ref, bestTrans, bestMatch);
		assertEquals(rmsd, rmsdCopy, 0.001);
		assertEquals(ref.getNumAtoms(), bestMatch.size());
	}

	@Test
	public void testBestRMSConjugatedGroups() {
		ROMol mol = RWMol.MolFromSmiles(
				"CCC(=O)[O-] |(-1.11,0.08,-0.29;0.08,-0.18,0.58;1.34,0.03,-0.16;1.74,1.22,-0.32;2.06,-1.04,-0.66)|");
		ROMol qry = RWMol.MolFromSmiles(
				"CCC([O-])=O |(-1.11,0.08,-0.29;0.08,-0.18,0.58;1.34,0.03,-0.16;1.74,1.22,-0.32;2.06,-1.04,-0.66)|");
		assertEquals(0.0, RDKFuncs.getBestRMS(qry, mol), 0.001);
		assertEquals(0.747,
				RDKFuncs.getBestRMS(qry, mol, -1, -1, new Match_Vect_Vect(), 1000000, false),
				0.001);
	}

	// equivalent of Python's Chem.MultiConfMolFromSDF()
	private ROMol multiConfMolFromSDF(String name) {
		SDMolSupplier suppl = new SDMolSupplier(alignDataFile(name));
		ROMol mol = new ROMol(suppl.next());
		while (!suppl.atEnd()) {
			ROMol m = suppl.next();
			mol.addConformer(new Conformer(m.getConformer()), true);
		}
		return mol;
	}

	@Test
	public void testBestAlignmentParams() {
		SDMolSupplier suppl = new SDMolSupplier(alignDataFile("probe_mol.sdf"), true, false);
		suppl.moveTo(1);
		ROMol prb = suppl.next();
		ROMol ref = suppl.next();

		// Hs are ignored by default
		BestAlignmentParams params = new BestAlignmentParams();
		assertTrue(params.getIgnoreHs());
		assertEquals(1.8100, RDKFuncs.getBestRMS(new ROMol(prb), ref, params), 0.001);

		params.setIgnoreHs(false);
		assertEquals(2.43449, RDKFuncs.getBestRMS(new ROMol(prb), ref, params), 0.001);

		Transform3D bestTrans = new Transform3D();
		Match_Vect bestMatch = new Match_Vect();
		double rmsd = RDKFuncs.getBestAlignmentTransform(new ROMol(prb), ref, bestTrans,
				bestMatch, params);
		assertEquals(2.43449, rmsd, 0.001);
		assertEquals(ref.getNumAtoms(), bestMatch.size());
	}

	@Test
	public void testGetAllConformerBestRMSToRef() {
		ROMol prbMol = multiConfMolFromSDF("butane_prb.sdf");
		SDMolSupplier refSuppl = new SDMolSupplier(alignDataFile("butane_ref.sdf"));
		ROMol refMol = refSuppl.next();
		BestAlignmentParams params = new BestAlignmentParams();
		double[] expected = {0.19474, 0.86739, 0.87102, 0.35358, 0.35395};
		Double_Vect rmsds = RDKFuncs.getAllConformerBestRMSToRef(prbMol, refMol, params);
		assertEquals(expected.length, rmsds.size());
		for (int i = 0; i < expected.length; ++i) {
			assertEquals(expected[i], rmsds.get(i), 0.0001);
		}

		refMol = multiConfMolFromSDF("butane_ref.sdf");
		double[] expectedMulti = {0.19474, 0.86739, 0.87102, 0.35358, 0.35395,
				0.82243, 0.16809, 0.16859, 0.54966, 0.56173};
		rmsds = RDKFuncs.getAllConformerBestRMSToRef(prbMol, refMol, params);
		assertEquals(expectedMulti.length, rmsds.size());
		for (int i = 0; i < expectedMulti.length; ++i) {
			assertEquals(expectedMulti[i], rmsds.get(i), 0.0001);
		}
	}

	@Test
	public void testGetAllConformerBestRMS() {
		ROMol mol = multiConfMolFromSDF("symmetric.confs.sdf");
		long nconfs = mol.getNumConformers();
		assertTrue(nconfs > 1);
		Double_Vect origVals = RDKFuncs.getAllConformerBestRMS(mol);
		assertEquals(nconfs * (nconfs - 1) / 2, origVals.size());

		Double_Vect newVals = RDKFuncs.getAllConformerBestRMS(mol, 4);
		assertEquals(origVals.size(), newVals.size());
		for (int i = 0; i < origVals.size(); ++i) {
			assertEquals(origVals.get(i), newVals.get(i), 0.0001);
		}

		assertEquals(origVals.get(0), RDKFuncs.getBestRMS(mol, mol, 0, 1), 0.0001);
	}

	private void checkO3AOnRefE2(O3A.AtomTypeScheme atomTypes, double expectedScore,
			double expectedRMSD) {
		SDMolSupplier suppl = new SDMolSupplier(alignDataFile("ref_e2.sdf"), true, false);
		suppl.moveTo(48);
		ROMol refMol = suppl.next();
		suppl.reset();
		double cumScore = 0.0;
		double cumMsd = 0.0;
		int nMols = 0;
		while (!suppl.atEnd()) {
			ROMol prbMol = suppl.next();
			O3A o3a = new O3A(prbMol, refMol, atomTypes);
			cumScore += o3a.score();
			double rmsd = o3a.align();
			cumMsd += rmsd * rmsd;
			++nMols;
			o3a.delete();
		}
		cumMsd /= nMols;
		assertEquals(expectedScore, cumScore, 1.0);
		assertEquals(expectedRMSD, Math.sqrt(cumMsd), 0.001);
	}

	@Test
	public void testO3AMMFF() {
		checkO3AOnRefE2(O3A.AtomTypeScheme.MMFF94, 6942, 0.345);
	}

	@Test
	public void testO3ACrippen() {
		checkO3AOnRefE2(O3A.AtomTypeScheme.CRIPPEN, 4918, 0.304);
	}

	@Test
	public void testO3AMatchesAndTransform() {
		SDMolSupplier suppl = new SDMolSupplier(alignDataFile("ref_e2.sdf"), true, false);
		ROMol refMol = suppl.next();
		ROMol prbMol = suppl.next();
		O3A o3a = new O3A(prbMol, refMol);
		Match_Vect matches = o3a.matches();
		assertTrue(matches.size() > 0);
		assertEquals(matches.size(), o3a.weights().size());
		Transform3D trans = new Transform3D();
		double transRmsd = o3a.trans(trans);
		assertEquals(transRmsd, o3a.align(), 1e-6);
	}

	@Test(expected = GenericRDKitException.class)
	public void testO3AMissingMMFFParams() {
		ROMol m1 = RWMol.MolFromSmiles("c1ccccc1Cl");
		DistanceGeom.EmbedMolecule(m1);
		ROMol m2 = RWMol.MolFromSmiles("c1ccccc1B(O)O");
		DistanceGeom.EmbedMolecule(m2);
		new O3A(m1, m2);
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.AlignTests");
	}

}
