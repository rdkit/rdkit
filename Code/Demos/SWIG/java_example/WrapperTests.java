// $Id$
// Copyright (C) 2008 Greg Landrum
//
// @@  All Rights Reserved @@
import org.junit.*;
import static org.junit.Assert.*;

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

    static {
        try {
            System.loadLibrary("RDKFuncs");
        } catch (UnsatisfiedLinkError e) {
            System.err.println("Native code library failed to load. See the chapter on Dynamic Linking Problems in the SWIG Java documentation for help.\n" + e);
            System.exit(1);
        }
    }

    public static void main(String args[]) {
      org.junit.runner.JUnitCore.main("WrapperTests");
    }
}


