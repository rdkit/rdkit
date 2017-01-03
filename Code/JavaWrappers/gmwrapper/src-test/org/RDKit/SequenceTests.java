/*
 *
 *  Copyright (c) 2016, Greg Landrum
 *
 *   @@ All Rights Reserved @@
 *  This file is part of the RDKit.
 *  The contents are covered by the terms of the BSD license
 *  which is included in the file license.txt, found at the root
 *  of the RDKit source tree.
*/

package org.RDKit;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.Test;

public class SequenceTests extends GraphMolTest {
	@Test
	public void testSequence1() {
      ROMol m = RWMol.MolFromSequence("CYIQNCPLG");
      AtomMonomerInfo mi=new AtomMonomerInfo(m.getAtomWithIdx(0).getMonomerInfo());
      assert(mi instanceof AtomPDBResidueInfo);
      String seq = new String(m.MolToSequence());
      assertEquals(seq,"CYIQNCPLG");
      String fasta = new String(m.MolToFASTA());
      assertEquals(fasta,">\nCYIQNCPLG\n");
      String helm = new String(m.MolToHELM());
      assertEquals(helm,"PEPTIDE1{C.Y.I.Q.N.C.P.L.G}$$$$");
	}
  @Test
	public void testSequence2() {
      ROMol m = RWMol.MolFromFASTA(">\nCYIQNCPLG\n");
      AtomMonomerInfo mi=new AtomMonomerInfo(m.getAtomWithIdx(0).getMonomerInfo());
      assert(mi instanceof AtomPDBResidueInfo);
      String seq = new String(m.MolToSequence());
      assertEquals(seq,"CYIQNCPLG");
      String fasta = new String(m.MolToFASTA());
      assertEquals(fasta,">\nCYIQNCPLG\n");
      String helm = new String(m.MolToHELM());
      assertEquals(helm,"PEPTIDE1{C.Y.I.Q.N.C.P.L.G}$$$$");
	}
  @Test
	public void testSequence3() {
      ROMol m = RWMol.MolFromHELM("PEPTIDE1{C.Y.I.Q.N.C.P.L.G}$$$$\n");
      AtomMonomerInfo mi=new AtomMonomerInfo(m.getAtomWithIdx(0).getMonomerInfo());
      assert(mi instanceof AtomPDBResidueInfo);
      String seq = new String(m.MolToSequence());
      assertEquals(seq,"CYIQNCPLG");
      String fasta = new String(m.MolToFASTA());
      assertEquals(fasta,">\nCYIQNCPLG\n");
      String helm = new String(m.MolToHELM());
      assertEquals(helm,"PEPTIDE1{C.Y.I.Q.N.C.P.L.G}$$$$");
	}
  @Test
	public void testSequence4() {
      ROMol m = RWMol.MolFromSequence("CGCGAATTACCGCG",false,6);
      AtomMonomerInfo mi=new AtomMonomerInfo(m.getAtomWithIdx(0).getMonomerInfo());
      assert(mi instanceof AtomPDBResidueInfo);
      String seq = new String(m.MolToSequence());
      assertEquals(seq,"CGCGAATTACCGCG");
      String fasta = new String(m.MolToFASTA());
      assertEquals(fasta,">\nCGCGAATTACCGCG\n");
      String helm = new String(m.MolToHELM());
      assertEquals(helm,"RNA1{[dR](C)P.[dR](G)P.[dR](C)P.[dR](G)P.[dR](A)P.[dR](A)P.[dR](T)P.[dR](T)P.[dR](A)P.[dR](C)P.[dR](C)P.[dR](G)P.[dR](C)P.[dR](G)}$$$$");
	}
  public void testSequence5() {
      ROMol m = RWMol.MolFromSequence("CGCGAAUUACCGCG",false,2);
      AtomMonomerInfo mi=new AtomMonomerInfo(m.getAtomWithIdx(0).getMonomerInfo());
      assert(mi instanceof AtomPDBResidueInfo);
      String seq = new String(m.MolToSequence());
      assertEquals(seq,"CGCGAAUUACCGCG");
      String fasta = new String(m.MolToFASTA());
      assertEquals(fasta,">\nCGCGAAUUACCGCG\n");
      String helm = new String(m.MolToHELM());
      assertEquals(helm,"RNA1{R(C)P.R(G)P.R(C)P.R(G)P.R(A)P.R(A)P.R(U)P.R(U)P.R(A)P.R(C)P.R(C)P.R(G)P.R(C)P.R(G)}$$$$");
	}


	public static void main(String args[]) {
            org.junit.runner.JUnitCore.main("org.RDKit.SequenceTests");
	}

}
