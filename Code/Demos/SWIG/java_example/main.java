// $Id$
// Copyright (C) 2008 Greg Landrum
//
// @@  All Rights Reserved @@

import org.RDKit.*;

public class main {
  static {
    try {
	System.loadLibrary("RDKFuncs");
    } catch (UnsatisfiedLinkError e) {
      System.err.println("Native code library failed to load. See the chapter on Dynamic Linking Problems in the SWIG Java documentation for help.\n" + e);
      System.exit(1);
    }
  }
  public static void main(String argv[]) {
      String smiles="c1ccccc1";
      ROMol mol= RDKFuncs.MolFromSmiles(smiles);
      long nAtoms = mol.getNumAtoms();
      System.out.println("nAtoms: " + nAtoms);
      System.out.println("smi: " + RDKFuncs.MolToSmiles(mol));
      System.out.println("atom: " + mol.getAtomWithIdx(0).getAtomicNum());
      System.out.println("bond: " + mol.getBondWithIdx(0).getIdx());
      System.out.println("bond: " + mol.getBondWithIdx(0).getBondType());
      System.out.println("hss: " + mol.hasSubstructMatch(RDKFuncs.MolFromSmarts("c")));
      System.out.println("hss: " + mol.hasSubstructMatch(RDKFuncs.MolFromSmarts("C")));
      System.out.println("hss: " + mol.hasSubstructMatch(RDKFuncs.MolFromSmarts("C")));
      RingInfo rI= mol.getRingInfo();
      System.out.println("ri: " + rI.isAtomInRingOfSize(0,6) + " " + rI.isAtomInRingOfSize(0,5) );
  }
}
