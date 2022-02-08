
//
// Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//
package org.RDKit;

import static org.junit.Assert.*;
import org.junit.Test;

public class MolEnumeratorTests extends  GraphMolTest {

    @Test
    public void linkNodes() {
        RWMol mol = RWMol.MolFromMolBlock("one linknode\n"+
"  Mrv2007 06222005102D          \n"+
"\n"+
"  0  0  0     0  0            999 V3000\n"+
"M  V30 BEGIN CTAB\n"+
"M  V30 COUNTS 6 6 0 0 0\n"+
"M  V30 BEGIN ATOM\n"+
"M  V30 1 C 8.25 12.1847 0 0\n"+
"M  V30 2 C 6.9164 12.9547 0 0\n"+
"M  V30 3 C 6.9164 14.4947 0 0\n"+
"M  V30 4 C 9.5836 14.4947 0 0\n"+
"M  V30 5 C 9.5836 12.9547 0 0\n"+
"M  V30 6 O 8.25 10.6447 0 0\n"+
"M  V30 END ATOM\n"+
"M  V30 BEGIN BOND\n"+
"M  V30 1 1 1 2\n"+
"M  V30 2 1 2 3\n"+
"M  V30 3 1 4 5\n"+
"M  V30 4 1 1 5\n"+
"M  V30 5 1 3 4\n"+
"M  V30 6 1 1 6\n"+
"M  V30 END BOND\n"+
"M  V30 LINKNODE 1 4 2 1 2 1 5\n"+
"M  V30 END CTAB\n"+
"M  END");
        MolEnumeratorParams ps = RDKFuncs.getLinkNodeParams();
        MolBundle bndl = RDKFuncs.enumerate(mol,ps);
        assertEquals(bndl.size(),4);
        bndl.delete();        
        mol.delete();
    }

    @Test
    public void positionVariation() {
        RWMol mol = RWMol.MolFromMolBlock("one variable attachment\n"+
"  Mrv2007 06232015292D          \n"+
"\n"+
"  0  0  0     0  0            999 V3000\n"+
"M  V30 BEGIN CTAB\n"+
"M  V30 COUNTS 9 8 0 0 0\n"+
"M  V30 BEGIN ATOM\n"+
"M  V30 1 C -1.7083 2.415 0 0\n"+
"M  V30 2 C -3.042 1.645 0 0\n"+
"M  V30 3 C -3.042 0.105 0 0\n"+
"M  V30 4 N -1.7083 -0.665 0 0\n"+
"M  V30 5 C -0.3747 0.105 0 0\n"+
"M  V30 6 C -0.3747 1.645 0 0\n"+
"M  V30 7 * -0.8192 1.3883 0 0\n"+
"M  V30 8 O -0.8192 3.6983 0 0\n"+
"M  V30 9 C 0.5145 4.4683 0 0\n"+
"M  V30 END ATOM\n"+
"M  V30 BEGIN BOND\n"+
"M  V30 1 1 1 2\n"+
"M  V30 2 2 2 3\n"+
"M  V30 3 1 3 4\n"+
"M  V30 4 2 4 5\n"+
"M  V30 5 1 5 6\n"+
"M  V30 6 2 1 6\n"+
"M  V30 7 1 7 8 ENDPTS=(3 1 5 6) ATTACH=ANY\n"+
"M  V30 8 1 8 9\n"+
"M  V30 END BOND\n"+
"M  V30 END CTAB\n"+
"M  END");
        MolEnumeratorParams ps = RDKFuncs.getPositionVariationParams();
        MolBundle bndl = RDKFuncs.enumerate(mol,ps);
        assertEquals(bndl.size(),3);
        
        mol.delete();
    }

    @Test
    public void enumerateMolecule() {
        RWMol mol = RWMol.MolFromMolBlock("one linknode\n"+
"  Mrv2007 06222005102D          \n"+
"\n"+
"  0  0  0     0  0            999 V3000\n"+
"M  V30 BEGIN CTAB\n"+
"M  V30 COUNTS 6 6 0 0 0\n"+
"M  V30 BEGIN ATOM\n"+
"M  V30 1 C 8.25 12.1847 0 0\n"+
"M  V30 2 C 6.9164 12.9547 0 0\n"+
"M  V30 3 C 6.9164 14.4947 0 0\n"+
"M  V30 4 C 9.5836 14.4947 0 0\n"+
"M  V30 5 C 9.5836 12.9547 0 0\n"+
"M  V30 6 O 8.25 10.6447 0 0\n"+
"M  V30 END ATOM\n"+
"M  V30 BEGIN BOND\n"+
"M  V30 1 1 1 2\n"+
"M  V30 2 1 2 3\n"+
"M  V30 3 1 4 5\n"+
"M  V30 4 1 1 5\n"+
"M  V30 5 1 3 4\n"+
"M  V30 6 1 1 6\n"+
"M  V30 END BOND\n"+
"M  V30 LINKNODE 1 4 2 1 2 1 5\n"+
"M  V30 END CTAB\n"+
"M  END");
        MolBundle bndl = RDKFuncs.enumerate(mol);
        assertEquals(bndl.size(),4);
        bndl.delete();        
        mol.delete();
    }

    @Test
    public void enumerateEmpty() {
        RWMol mol = RWMol.MolFromMolBlock("no linknodes\n"+
"  Mrv2007 06222005102D          \n"+
"\n"+
"  0  0  0     0  0            999 V3000\n"+
"M  V30 BEGIN CTAB\n"+
"M  V30 COUNTS 6 6 0 0 0\n"+
"M  V30 BEGIN ATOM\n"+
"M  V30 1 C 8.25 12.1847 0 0\n"+
"M  V30 2 C 6.9164 12.9547 0 0\n"+
"M  V30 3 C 6.9164 14.4947 0 0\n"+
"M  V30 4 C 9.5836 14.4947 0 0\n"+
"M  V30 5 C 9.5836 12.9547 0 0\n"+
"M  V30 6 O 8.25 10.6447 0 0\n"+
"M  V30 END ATOM\n"+
"M  V30 BEGIN BOND\n"+
"M  V30 1 1 1 2\n"+
"M  V30 2 1 2 3\n"+
"M  V30 3 1 4 5\n"+
"M  V30 4 1 1 5\n"+
"M  V30 5 1 3 4\n"+
"M  V30 6 1 1 6\n"+
"M  V30 END BOND\n"+
"M  V30 END CTAB\n"+
"M  END");
        MolBundle bndl = RDKFuncs.enumerate(mol);
        assertEquals(bndl.size(),0);
        bndl.delete();        
        mol.delete();
    }


    public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("org.RDKit.MolEnumeratorTests");
    }

}
