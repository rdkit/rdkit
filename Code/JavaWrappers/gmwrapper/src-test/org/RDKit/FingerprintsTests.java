/*
 * $Id: FingerprintsTests.java 131 2011-01-20 22:01:29Z ebakke $
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

public class FingerprintsTests extends GraphMolTest {

	public void compareVectors(Int_Vect v1, Int_Vect v2) {
		assertEquals(v1.size(), v2.size());
		for (int i = 0; i < v1.size(); i++) {
			assertEquals(v1.get(i),v2.get(i));
		}
	}

	private void checkTargets(String smi, String[] matches) {
		ROMol m = RWMol.MolFromSmiles(smi);
		ExplicitBitVect fp1 = RDKFuncs.RDKFingerprintMol(m, 2, 7, 9192, 4, false);
		for (String match : matches) {
			ROMol m2 = RWMol.MolFromSmiles(match);
			ExplicitBitVect fp2 = RDKFuncs.RDKFingerprintMol(m2, 2, 7, 9192, 4, false);
			Double_Vect v = RDKFuncs.OnBitProjSimilarity(fp2, fp1);
			assertEquals("substuct " + match + " is not properly contained in " + smi,1.000,v.get(0),defaultDoubleTol);
		}
	}

	// check containing mols, no Hs, no valence
	@Test
	public void test4 () {
		String smi = "CCC(O)C(=O)O";
		String[] matches = new String[] {"CCC","OCC","OCC=O","OCCO","CCCC","OC=O","CC(O)C"};
		checkTargets(smi,matches);
	}

	// check containing mols, use Hs, no valence
	@Test
	public void test5() {
		String smi = "CCC(O)C(=O)O";
		String[] matches = new String[] {"O[CH-][CH2-]","O[CH-][C-]=O"};
		checkTargets(smi,matches);
	}

	// check that the bits in a signature of size N which has been folded in half
    // are the same as those in a signature of size N/2
	@Test
	public void test6() {
		String[] smis = new String[] {
				"CCC(O)C(=O)O","c1ccccc1","C1CCCCC1","C1NCCCC1","CNCNCNC"
		};
		for (String smi : smis) {
			ROMol m = RWMol.MolFromSmiles(smi);
			ExplicitBitVect fp1 = RDKFuncs.RDKFingerprintMol(m, 2, 7, 4096);
			ExplicitBitVect fp3 = RDKFuncs.RDKFingerprintMol(m, 2, 7, 2048);
			ExplicitBitVect fp2 = RDKFuncs.FoldFingerprint(fp1, 2);
			compareVectors(fp2.getOnBits(), fp3.getOnBits());
			fp2 = RDKFuncs.FoldFingerprint(fp2, 2);
			fp3 = RDKFuncs.RDKFingerprintMol(m, 2, 7, 1024);
			compareVectors(fp2.getOnBits(), fp3.getOnBits());
			fp2 = RDKFuncs.FoldFingerprint(fp1, 4);
			compareVectors(fp2.getOnBits(), fp3.getOnBits());
		}
	}

	@Test
	public void test7() {
		String smi1 = "c1ccccc1";
                String smi2 = "c1ccccn1";
		ROMol m1 = RWMol.MolFromSmiles(smi1);
                ROMol m2 = RWMol.MolFromSmiles(smi2);
                ExplicitBitVect fp1 = RDKFuncs.MACCSFingerprintMol(m1);
                ExplicitBitVect fp2 = RDKFuncs.MACCSFingerprintMol(m2);
                assertEquals(RDKFuncs.DiceSimilarity(fp1,fp2),0.5454,0.001);
	}

	@Test
	public void test8() {
		ReactionFingerprintParams params = new ReactionFingerprintParams();
		params.setFpType(FingerprintType.PatternFP);
		params.setFpSize(4096);
		{
			String smi1 = "C1CCCCC1>>C1CCNCC1";
	    String smi2 = "C1CCCCC1>>C1CCNCC1";
			ChemicalReaction r1 = ChemicalReaction.ReactionFromSmarts(smi1,true);
			ChemicalReaction r2 = ChemicalReaction.ReactionFromSmarts(smi2,true);
	    ExplicitBitVect fp1 = RDKFuncs.StructuralFingerprintChemReaction(r1,params);
	    ExplicitBitVect fp2 = RDKFuncs.StructuralFingerprintChemReaction(r2,params);
	    assertTrue(RDKFuncs.AllProbeBitsMatch(fp1,fp2));
	}
	{
		String smi1 = "C1CCCCC1>>C1CCNCC1";
		String smi2 = "C1CCCCC1>>C1CCOCC1";
		ChemicalReaction r1 = ChemicalReaction.ReactionFromSmarts(smi1,true);
		ChemicalReaction r2 = ChemicalReaction.ReactionFromSmarts(smi2,true);
		ExplicitBitVect fp1 = RDKFuncs.StructuralFingerprintChemReaction(r1,params);
		ExplicitBitVect fp2 = RDKFuncs.StructuralFingerprintChemReaction(r2,params);
		assertFalse(RDKFuncs.AllProbeBitsMatch(fp1,fp2));
	}
	{
		String smi1 = "C1CCCCC1>>C1CCNCC1";
		String smi2 = ">>C1CCNCC1";
		ChemicalReaction r1 = ChemicalReaction.ReactionFromSmarts(smi1,true);
		ChemicalReaction r2 = ChemicalReaction.ReactionFromSmarts(smi2,true);
		ExplicitBitVect fp1 = RDKFuncs.StructuralFingerprintChemReaction(r1,params);
		ExplicitBitVect fp2 = RDKFuncs.StructuralFingerprintChemReaction(r2,params);
		assertFalse(RDKFuncs.AllProbeBitsMatch(fp1,fp2));
		assertTrue(RDKFuncs.AllProbeBitsMatch(fp2,fp1));
	}
}

    @Test
    public void testToByteArray() {
        String smiles = "Cc2nc1ccccc1o2";
        ROMol mol = RWMol.MolFromSmiles(smiles);
        ExplicitBitVect fp1 = RDKFuncs.PatternFingerprintMol(mol, 2048);
        byte[] fpBytes = fp1.toByteArray();
        ExplicitBitVect fp2 = ExplicitBitVect.fromByteArray(fpBytes);
        Int_Vect fp1Bits = fp1.getOnBits();
        Int_Vect fp2Bits = fp2.getOnBits();
        compareVectors(fp1Bits, fp2Bits);
    }

    @Test
    public void testPatternFpsForBundles() {
		ROMol q1 = RWMol.MolFromSmiles("OCCO");
		ROMol q2 = RWMol.MolFromSmiles("OCCCO");
        ExplicitBitVect qfp1 = RDKFuncs.PatternFingerprintMol(q1, 2048);
        ExplicitBitVect qfp2 = RDKFuncs.PatternFingerprintMol(q2, 2048);
		MolBundle bndl = new MolBundle();
		bndl.addMol(q1);
		bndl.addMol(q2);
		
		ExplicitBitVect bndlfp = RDKFuncs.PatternFingerprintMol(bndl,2048);
        Int_Vect fpBits = bndlfp.getOnBits();
		assertTrue(fpBits.size()!=0);
		for(int i=0;i<fpBits.size();i++){
			int idx = fpBits.get(i);
			assertTrue(qfp1.getBit(idx));
			assertTrue(qfp2.getBit(idx));
		}
		
    }

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.FingerprintsTests");
	}
}
