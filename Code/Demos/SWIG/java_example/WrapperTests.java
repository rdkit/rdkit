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
    }
    @Test public void testBonds() {
        assertTrue(mol1.getBondWithIdx(0).getBondType()==Bond.BondType.AROMATIC);
    }
    @Test public void testSmilesWrite() {
        String smi=RDKFuncs.MolToSmiles(mol1);
        assertEquals(smi,"c1ccccc1",smi);
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


