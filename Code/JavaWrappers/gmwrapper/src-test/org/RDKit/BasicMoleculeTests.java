/*
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

import org.junit.*;

public class BasicMoleculeTests extends GraphMolTest {

	private ROMol mol1;
	@Before public void setUp() {
		String smiles="c1ccccc1";
		mol1 = RWMol.MolFromSmiles(smiles);
	}
	@Test public void testBasics() {
		assertEquals(6,mol1.getNumAtoms());
		assertEquals(6,mol1.getNumBonds());
	}
	@Test public void testSmilesWrite() {
		String smi = mol1.MolToSmiles();
		assertEquals(smi,smi,"c1ccccc1");
		String smiRDKFuncs = RDKFuncs.MolToSmiles(mol1);
		assertEquals(smiRDKFuncs,smiRDKFuncs,smi);
	}
	@Test public void testAtoms() {
		assertEquals( mol1.getAtomWithIdx(0).getAtomicNum(),6);
		assertFalse(mol1.hasAtomBookmark(1));
		mol1.setAtomBookmark(mol1.getAtomWithIdx(0),1);
		assertTrue(mol1.hasAtomBookmark(1));
	}
	@Test public void testBonds() {
		assertEquals(Bond.BondType.AROMATIC,mol1.getBondWithIdx(0).getBondType());
	}
	@Test public void testReactionBasics() {
		ChemicalReaction rxn;
		rxn=ChemicalReaction.ReactionFromSmarts("[OH][C:1]=[O:2].[N!H0:3]>>[N:3][C:1]=[O:2]");
		assertEquals(2,rxn.getNumReactantTemplates());
		assertEquals(1,rxn.getNumProductTemplates());
		ROMol r1,r2;
		r1=RWMol.MolFromSmiles("CC(=O)O");
		r2=RWMol.MolFromSmiles("ClCN");
		assertEquals(4,r1.getNumAtoms());
		assertEquals(3,r2.getNumAtoms());
		ROMol_Vect rs= new ROMol_Vect(2);
		rs.set(0,r1);
		rs.set(1,r2);
		ROMol_Vect_Vect ps;
		ps=rxn.runReactants(rs);

		assertFalse(ps.isEmpty());
		assertEquals(1,ps.size());
		assertFalse(ps.get(0).isEmpty());
		assertEquals(1,ps.get(0).size());

		assertEquals(4,r1.getNumAtoms());
		assertEquals(3,r2.getNumAtoms());

		assertEquals(6,ps.get(0).get(0).getNumAtoms());
	}

	@Test public void testSubstruct1() {
		ROMol p;
		Match_Vect mv;
		p = RWMol.MolFromSmarts("c");
		assertTrue(mol1.hasSubstructMatch(p));
		mv=mol1.getSubstructMatch(p);
		assertEquals(1,mv.size());
		assertEquals(0,mv.get(0).getFirst());
		assertEquals(0,mv.get(0).getSecond());
	}
	@Test public void testSubstruct2() {
		ROMol p;
		Match_Vect mv;
		p = RWMol.MolFromSmarts("C");
		assertFalse(mol1.hasSubstructMatch(p));
		mv=mol1.getSubstructMatch(p);
		assertEquals(0,mv.size());
	}
	@Test public void testSubstruct3() {
		ROMol p;
		Match_Vect mv;
		ROMol m2;
		m2 = RWMol.MolFromSmiles("NC(=O)CC");
		p = RWMol.MolFromSmarts("CN");
		mv=m2.getSubstructMatch(p);
		assertEquals(2,mv.size());
		assertEquals(0,mv.get(0).getFirst());
		assertEquals(1,mv.get(0).getSecond());
		assertEquals(1,mv.get(1).getFirst());
		assertEquals(0,mv.get(1).getSecond());
	}
	@Test public void testSubstruct4() {
		ROMol p;
		Match_Vect_Vect mvv;
		ROMol m2;
		m2 = RWMol.MolFromSmiles("NC(=O)CC");
		p = RWMol.MolFromSmarts("CN");
		mvv=m2.getSubstructMatches(p);
		assertEquals(1,mvv.size());
		assertEquals(2,mvv.get(0).size());
		assertEquals(0,mvv.get(0).get(0).getFirst());
		assertEquals(1,mvv.get(0).get(0).getSecond());
		assertEquals(1,mvv.get(0).get(1).getFirst());
		assertEquals(0,mvv.get(0).get(1).getSecond());
	}
	@Test public void testSubstruct5() {
		ROMol p;
		Match_Vect_Vect mvv;
		ROMol m2;
		m2 = RWMol.MolFromSmiles("NC(=O)NCC");
		p = RWMol.MolFromSmarts("[$(C=O)]N");
		mvv=m2.getSubstructMatches(p);
		assertEquals(2,mvv.size());
		assertEquals(2,mvv.get(0).size());
		assertEquals(0,mvv.get(0).get(0).getFirst());
		assertEquals(1,mvv.get(0).get(0).getSecond());
		assertEquals(1,mvv.get(0).get(1).getFirst());
		assertEquals(0,mvv.get(0).get(1).getSecond());
		assertEquals(2,mvv.get(1).size());
		assertEquals(0,mvv.get(1).get(0).getFirst());
		assertEquals(1,mvv.get(1).get(0).getSecond());
		assertEquals(1,mvv.get(1).get(1).getFirst());
		assertEquals(3,mvv.get(1).get(1).getSecond());
	}
	@Test public void testSubstructParams1() {
		ROMol p;
		Match_Vect_Vect mvv;
		ROMol m2;
		m2 = RWMol.MolFromSmiles("C[C@](F)(Cl)Br");
		p = RWMol.MolFromSmiles("C[C@@](F)(Cl)Br");
		mvv=m2.getSubstructMatches(p);
		assertEquals(1,mvv.size());
		assertEquals(5,mvv.get(0).size());
		SubstructMatchParameters params = new SubstructMatchParameters();
		params.setUseChirality(true);
		mvv=m2.getSubstructMatches(p,params);
		assertEquals(0,mvv.size());
	}

	@Test public void testFingerprints1() {
		ROMol m1,m2;
		m1 = RWMol.MolFromSmiles("C1=CC=CC=C1");
		m2 = RWMol.MolFromSmiles("C1=CC=CC=N1");
		ExplicitBitVect fp1,fp2;
		fp1=RDKFuncs.RDKFingerprintMol(m1);
		fp2=RDKFuncs.RDKFingerprintMol(m1);
		assertEquals( 0.0,RDKFuncs.TanimotoSimilarity(fp1,fp2),1.0);
		fp2=RDKFuncs.RDKFingerprintMol(m2);
		assertTrue(RDKFuncs.TanimotoSimilarity(fp1,fp2)<1.0);
		assertTrue(RDKFuncs.TanimotoSimilarity(fp1,fp2)>0.0);
	}
	@Test public void testFingerprints2() {
		ROMol m1,m2;
		m1 = RWMol.MolFromSmiles("C1=CC=CC=C1");
		m2 = RWMol.MolFromSmiles("C1=CC=CC=N1");
		SparseIntVectu32 fp1,fp2;
		fp1=RDKFuncs.MorganFingerprintMol(m1,2);
		fp2=RDKFuncs.MorganFingerprintMol(m1,2);
		assertEquals(0.0,RDKFuncs.DiceSimilarity(fp1,fp2),1.0);
		fp2=RDKFuncs.MorganFingerprintMol(m2,2);
		assertTrue(RDKFuncs.DiceSimilarity(fp1,fp2)<1.0);
		assertTrue(RDKFuncs.DiceSimilarity(fp1,fp2)>0.0);
		UInt_Pair_Vect v1=fp1.getNonzero();
		assertTrue(v1.size()>0);
		UInt_Pair_Vect v2=fp2.getNonzero();
		assertTrue(v2.size()>0);
		assertTrue(v2.size()>v1.size());
	}
	@Test public void testErrorHandling_BadParse() {
		ROMol m1 = RWMol.MolFromSmiles("C1CC");
		assertEquals(null,m1);
	}
	@Test(expected=MolSanitizeException.class)
	public void testErrorHandling_BadKekulize() {
		RWMol.MolFromSmiles("c1cc1");
	}
	@Test(expected=MolSanitizeException.class)
	public void testErrorHandling_BadKekulize2() {
            RWMol mol=RWMol.MolFromSmiles("c1cc1",0,false);
            RDKFuncs.sanitizeMol(mol);
	}
	@Test public void testErrorHandling_BadKekulize3() {
            RWMol mol=RWMol.MolFromSmiles("c1cc1",0,false);
            int ops=SanitizeFlags.SANITIZE_ALL.swigValue();
            int failedAt=RDKFuncs.sanitizeMol(mol,ops);
            assertEquals(failedAt,SanitizeFlags.SANITIZE_KEKULIZE.swigValue());
            ops ^= SanitizeFlags.SANITIZE_KEKULIZE.swigValue();
            failedAt=RDKFuncs.sanitizeMol(mol,ops);
            assertEquals(failedAt,0);
	}

	@Test public void testAddAtomsAndBonds() {
		RWMol mol = new RWMol();
		assertEquals(mol.getNumBonds(),0);
		Atom carbon = new Atom(6);
		mol.addAtom(carbon,false);
		mol.addAtom(carbon,false);
		assertEquals(mol.getNumAtoms(),2);
		Bond single = new Bond(Bond.BondType.SINGLE);
		single.setOwningMol(mol);
		single.setBeginAtomIdx(1);
		single.setEndAtomIdx(0);
		mol.addBond(single);
		assertEquals(mol.getNumBonds(),1);
		assertEquals(mol.getBondWithIdx(0).getBeginAtomIdx(),1);
		assertEquals(mol.getBondWithIdx(0).getEndAtomIdx(),0);

		// the molecule has copied the bond, so if we change it, there's
		// no impact on the molecule itself:
		single.setBeginAtomIdx(0);
		single.setEndAtomIdx(1);
		assertEquals(mol.getBondWithIdx(0).getBeginAtomIdx(),1);
		assertEquals(mol.getBondWithIdx(0).getEndAtomIdx(),0);

	}

	@Test public void testGetAtomNeighbors() {
		String smiles="CC(C)C";
		ROMol mol = RWMol.MolFromSmiles(smiles);
		assertEquals(mol.getAtomNeighbors(mol.getAtomWithIdx(0)).size(),1);
		assertEquals(mol.getAtomNeighbors(mol.getAtomWithIdx(1)).size(),3);
		assertEquals(mol.getAtomNeighbors(mol.getAtomWithIdx(2)).size(),1);
		assertEquals(mol.getAtomBonds(mol.getAtomWithIdx(0)).size(),1);
		assertEquals(mol.getAtomBonds(mol.getAtomWithIdx(1)).size(),3);
		assertEquals(mol.getAtomBonds(mol.getAtomWithIdx(2)).size(),1);
	}

/*	@Test -- the contents of this test now in UnitTestPickling, testIssue219
	public void testConformer(){
	}
*/
	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.BasicMoleculeTests");
	}
}
