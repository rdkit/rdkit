
//
// Created by Gareth Jones on 6/1/2020.
//
// Copyright 2020 Schrodinger, Inc
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.


package org.RDKit;

import static org.junit.Assert.*;
import org.junit.Test;

public class TautomerQueryTests extends  GraphMolTest {

    @Test
    public void basicOperations() {
        RWMol mol = RWMol.MolFromSmiles("O=C1CCCCC1");
        TautomerQuery tautomerQuery = TautomerQuery.fromMol(mol);
        assertEquals(2, tautomerQuery.getTautomers().size());
        Sizet_Vect modifiedAtoms = tautomerQuery.getModifiedAtoms();
        assertEquals(3, modifiedAtoms.size());


        RWMol target = RWMol.MolFromSmiles("OC1=CCCC(CC)C1");
        boolean match = tautomerQuery.isSubstructOf(target);
        assertTrue(match);
        Match_Vect_Vect allMatches = tautomerQuery.substructOf(target);
        assertEquals(1, allMatches.size());
        Match_Vect atomMatches = allMatches.get(0);
        assertEquals(7, atomMatches.size());

        SubstructMatchParameters params = new SubstructMatchParameters();
        ROMol_Vect matchingTautomers = new ROMol_Vect();
        allMatches = tautomerQuery.substructOf(target, params, matchingTautomers);
        assertEquals(1, allMatches.size());
        assertEquals(1, matchingTautomers.size());

        match = target.hasSubstructMatch(matchingTautomers.get(0));
        assertTrue(match);

        ExplicitBitVect tautomerFingerprint = tautomerQuery.patternFingerprintTemplate();
        ExplicitBitVect targetFingerprint = TautomerQuery.patternFingerprintTarget(target);
        boolean matchingFingerprints =  RDKFuncs.AllProbeBitsMatch(tautomerFingerprint, targetFingerprint);
        assertTrue(matchingFingerprints);

        tautomerFingerprint.delete();
        targetFingerprint.delete();
    }

    public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("org.RDKit.TautomerQueryTests");
    }

}
