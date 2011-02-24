// $Id$
// Copyright (C) 2008-2010 Greg Landrum
//
// @@  All Rights Reserved @@
//package org.RDKit;

import org.junit.*;
import static org.junit.Assert.*;
import org.RDKit.*;


public class WrapperTests {
    private ROMol mol1;
    @Before public void setUp() {
      String smiles="c1ccccc1";
      mol1 = RDKFuncs.MolFromSmiles(smiles);
    }
    @Test public void testBasics() {
        assertTrue(mol1.getNumAtoms()==6);
        assertTrue(mol1.getNumBonds()==6);
    }
    @Test public void testAtoms() {
        assertTrue(mol1.getAtomWithIdx(0).getAtomicNum()==6);
	assertFalse(mol1.hasAtomBookmark(1));
	mol1.setAtomBookmark(mol1.getAtomWithIdx(0),1);
	assertTrue(mol1.hasAtomBookmark(1));
	
    }
    @Test public void testBonds() {
        assertTrue(mol1.getBondWithIdx(0).getBondType()==Bond.BondType.AROMATIC);
    }
    @Test public void testSmilesWrite() {
        String smi=RDKFuncs.MolToSmiles(mol1);
        assertEquals(smi,"c1ccccc1",smi);
    }

    @Test public void testReactionBasics() {
	ChemicalReaction rxn;
	rxn=RDKFuncs.ReactionFromSmarts("[OH][C:1]=[O:2].[N!H0:3]>>[N:3][C:1]=[O:2]");
	assertTrue(rxn.getNumReactantTemplates()==2);
	assertTrue(rxn.getNumProductTemplates()==1);
	ROMol r1,r2;
	r1=RDKFuncs.MolFromSmiles("CC(=O)O");
	r2=RDKFuncs.MolFromSmiles("ClCN");
	assertTrue(r1.getNumAtoms()==4);
	assertTrue(r2.getNumAtoms()==3);
	ROMol_Vect rs= new ROMol_Vect(2);
        rs.set(0,r1);
        rs.set(1,r2);
	ROMol_Vect_Vect ps;
	ps=rxn.runReactants(rs);
        
	assertFalse(ps.isEmpty());
	assertTrue(ps.size()==1);
	assertFalse(ps.get(0).isEmpty());
	assertTrue(ps.get(0).size()==1);
        
	assertTrue(r1.getNumAtoms()==4);
	assertTrue(r2.getNumAtoms()==3);

	assertTrue(ps.get(0).get(0).getNumAtoms()==6);
        }
    @Test public void testSubstruct1() {
	ROMol p;
	Match_Vect mv;
	p = RDKFuncs.MolFromSmarts("c");
        assertTrue(mol1.hasSubstructMatch(p));
	mv=mol1.getSubstructMatch(p);
	assertTrue(mv.size()==1);
	assertTrue(mv.get(0).getFirst()==0);
	assertTrue(mv.get(0).getSecond()==0);
    }
    @Test public void testSubstruct2() {
	ROMol p;
	Match_Vect mv;
	p = RDKFuncs.MolFromSmarts("C");
        assertFalse(mol1.hasSubstructMatch(p));
	mv=mol1.getSubstructMatch(p);
	assertTrue(mv.size()==0);
    }
    @Test public void testSubstruct3() {
	ROMol p;
	Match_Vect mv;
	ROMol m2;
	m2 = RDKFuncs.MolFromSmiles("NC(=O)CC");
	p = RDKFuncs.MolFromSmarts("CN");
	mv=m2.getSubstructMatch(p);
	assertTrue(mv.size()==2);
	assertTrue(mv.get(0).getFirst()==0);
	assertTrue(mv.get(0).getSecond()==1);
	assertTrue(mv.get(1).getFirst()==1);
	assertTrue(mv.get(1).getSecond()==0);
    }	
    @Test public void testSubstruct4() {
	ROMol p;
	Match_Vect_Vect mvv;
	ROMol m2;
	m2 = RDKFuncs.MolFromSmiles("NC(=O)CC");
	p = RDKFuncs.MolFromSmarts("CN");
	mvv=m2.getSubstructMatches(p);
	assertTrue(mvv.size()==1);
	assertTrue(mvv.get(0).size()==2);
	assertTrue(mvv.get(0).get(0).getFirst()==0);
	assertTrue(mvv.get(0).get(0).getSecond()==1);
	assertTrue(mvv.get(0).get(1).getFirst()==1);
	assertTrue(mvv.get(0).get(1).getSecond()==0);
    }
    @Test public void testSubstruct5() {
	ROMol p;
	Match_Vect_Vect mvv;
	ROMol m2;
	m2 = RDKFuncs.MolFromSmiles("NC(=O)NCC");
	p = RDKFuncs.MolFromSmarts("[$(C=O)]N");
	mvv=m2.getSubstructMatches(p);
	assertTrue(mvv.size()==2);
	assertTrue(mvv.get(0).size()==2);
	assertTrue(mvv.get(0).get(0).getFirst()==0);
	assertTrue(mvv.get(0).get(0).getSecond()==1);
	assertTrue(mvv.get(0).get(1).getFirst()==1);
	assertTrue(mvv.get(0).get(1).getSecond()==0);
	assertTrue(mvv.get(1).size()==2);
	assertTrue(mvv.get(1).get(0).getFirst()==0);
	assertTrue(mvv.get(1).get(0).getSecond()==1);
	assertTrue(mvv.get(1).get(1).getFirst()==1);
	assertTrue(mvv.get(1).get(1).getSecond()==3);
    }
    @Test public void testFingerprints1() {
	ROMol m1,m2;
	m1 = RDKFuncs.MolFromSmiles("C1=CC=CC=C1");
	m2 = RDKFuncs.MolFromSmiles("C1=CC=CC=N1");
	ExplicitBitVect fp1,fp2;
	fp1=RDKFuncs.RDKFingerprintMol(m1);
	fp2=RDKFuncs.RDKFingerprintMol(m1);
	assertTrue(RDKFuncs.TanimotoSimilarityEBV(fp1,fp2)==1.0);
	fp2=RDKFuncs.RDKFingerprintMol(m2);
	assertTrue(RDKFuncs.TanimotoSimilarityEBV(fp1,fp2)<1.0);
	assertTrue(RDKFuncs.TanimotoSimilarityEBV(fp1,fp2)>0.0);
    }
    @Test public void testFingerprints2() {
	ROMol m1,m2;
	m1 = RDKFuncs.MolFromSmiles("C1=CC=CC=C1");
	m2 = RDKFuncs.MolFromSmiles("C1=CC=CC=N1");
	SparseIntVectu32 fp1,fp2;
	fp1=RDKFuncs.getMorganFingerprint(m1,2);
	fp2=RDKFuncs.getMorganFingerprint(m1,2);
	assertTrue(RDKFuncs.DiceSimilaritySIVu32(fp1,fp2)==1.0);
	fp2=RDKFuncs.getMorganFingerprint(m2,2);
	assertTrue(RDKFuncs.DiceSimilaritySIVu32(fp1,fp2)<1.0);
	assertTrue(RDKFuncs.DiceSimilaritySIVu32(fp1,fp2)>0.0);
	UInt_Pair_Vect v1=fp1.getNonzero();
	assertTrue(v1.size()>0);
	UInt_Pair_Vect v2=fp2.getNonzero();
	assertTrue(v2.size()>0);
	assertTrue(v2.size()>v1.size());
    }
    @Test public void testFingerprints3() {
	ROMol m1,m2;
	m1 = RDKFuncs.MolFromSmiles("C1=CC=CC=C1");
	m2 = RDKFuncs.MolFromSmiles("C1=CC=CC=N1");
	SparseIntVecti32 fp1,fp2;
	fp1=RDKFuncs.getAtomPairFingerprint(m1,2,6);
	fp2=RDKFuncs.getAtomPairFingerprint(m1,2,6);
	assertTrue(RDKFuncs.DiceSimilaritySIVi32(fp1,fp2)==1.0);
	fp2=RDKFuncs.getAtomPairFingerprint(m2,2,6);
        assertEquals(RDKFuncs.DiceSimilaritySIVi32(fp1,fp2),0.66667,.0001);
    }
    @Test public void testFingerprints4() {
	ROMol m1,m2;
	m1 = RDKFuncs.MolFromSmiles("C1=CC=CC=C1");
	m2 = RDKFuncs.MolFromSmiles("C1=CC=CC=N1");
	SparseIntVecti64 fp1,fp2;
	fp1=RDKFuncs.getTopologicalTorsionFingerprint(m1);
	fp2=RDKFuncs.getTopologicalTorsionFingerprint(m1);
	assertTrue(RDKFuncs.DiceSimilaritySIVi64(fp1,fp2)==1.0);
	fp2=RDKFuncs.getTopologicalTorsionFingerprint(m2);
        assertEquals(RDKFuncs.DiceSimilaritySIVi64(fp1,fp2),0.3333,.0001);
    }
    @Test public void testErrorHandling() {
	ROMol m1;
	m1 = RDKFuncs.MolFromSmiles("C1CC");
	assertTrue(m1==null);
	m1 = RDKFuncs.MolFromSmiles("c1cc1");
	assertTrue(m1==null);
	System.err.println("ok!");
        ChemicalReaction rxn=RDKFuncs.ReactionFromSmarts("OH][C:1]=[O:2].[N!H0:3]>>[N:3][C:1]=[O:2]");
        assertTrue(rxn==null);
    }
    @Test public void testPickling() {
	Int_Vect pkl=RDKFuncs.MolToBinary(mol1);
	ROMol m1 = RDKFuncs.MolFromBinary(pkl);
        assertTrue(m1.getNumAtoms()==6);
        assertTrue(m1.getNumBonds()==6);
    }
     @Test public void testReactionPickling() {
         ChemicalReaction rxn;
         rxn=RDKFuncs.ReactionFromSmarts("[OH][C:1]=[O:2].[N!H0:3]>>[N:3][C:1]=[O:2]");
         assertTrue(rxn.getNumReactantTemplates()==2);
         assertTrue(rxn.getNumProductTemplates()==1);
         Int_Vect pkl=RDKFuncs.RxnToBinary(rxn);
         ChemicalReaction rxn2=RDKFuncs.RxnFromBinary(pkl);
         assertTrue(rxn2.getNumReactantTemplates()==2);
         assertTrue(rxn2.getNumProductTemplates()==1);
    }
    @Test public void testReactionToSmarts() {
         ChemicalReaction rxn;
         rxn=RDKFuncs.ReactionFromSmarts("[OH][C:1]=[O:2].[N!H0:3]>>[N:3][C:1]=[O:2]");
         assertTrue(rxn.getNumReactantTemplates()==2);
         assertTrue(rxn.getNumProductTemplates()==1);
         String sma=RDKFuncs.ReactionToSmarts(rxn);
         ChemicalReaction rxn2;
         rxn2=RDKFuncs.ReactionFromSmarts(sma);
         assertTrue(rxn2.getNumReactantTemplates()==2);
         assertTrue(rxn2.getNumProductTemplates()==1);
    }


    /*@Test*/ public void testMemory() {
	for(int i=0;i<1000000;i++){
	    ROMol m1;
	    m1 = RDKFuncs.MolFromSmiles("C1=CC=CC=C1");
	    if((i%1000)==0){
		System.err.println("done: "+i);
	    }
	    m1.delete();
	}
    }
    /*@Test*/ public void testMemory2() {
	ChemicalReaction rxn;
	rxn=RDKFuncs.ReactionFromSmarts("[OH][C:1]=[O:2].[N!H0:3]>>[N:3][C:1]=[O:2]");
	assertTrue(rxn.getNumReactantTemplates()==2);
	assertTrue(rxn.getNumProductTemplates()==1);
	ROMol r1,r2;
	r1=RDKFuncs.MolFromSmiles("CC(=O)O");
	r2=RDKFuncs.MolFromSmiles("ClCN");
	assertTrue(r1.getNumAtoms()==4);
	assertTrue(r2.getNumAtoms()==3);
	for(int i=0;i<1000000;i++){
            ROMol_Vect rs= new ROMol_Vect(2);
            rs.set(0,r1);
            rs.set(1,r2);
            ROMol_Vect_Vect ps;
            ps=rxn.runReactants(rs);
	    if((i%1000)==0){
		System.err.println("done: "+i);
	    }
            ps.delete();
        }
    }
    /*@Test*/ public void testMemory3() {
        ROMol m1;
        m1 = RDKFuncs.MolFromSmiles("C1=CC=CC=C1C2CCC(C(=O)OCCC)CC2");
	for(int i=0;i<1000000;i++){
            /*SparseIntVecti64 fp1;
              fp1=RDKFuncs.getTopologicalTorsionFingerprint(m1);*/
            SparseIntVectu32 fp1;
            fp1=RDKFuncs.getMorganFingerprint(m1,2);
	    if((i%1000)==0){
		System.err.println("done: "+i);
	    }
            fp1.delete();
	}
    }
    static {
        try {
            System.loadLibrary("RDKFuncs");
        } catch (UnsatisfiedLinkError e) {
            System.err.println("Native code library failed to load. Make sure that libRDKFuncs.so is somewhere in your LD_LIBRARY_PATH.\n" + e);
            System.exit(1);
        }
    }

    public static void main(String args[]) {
      org.junit.runner.JUnitCore.main("WrapperTests");
    }
}


