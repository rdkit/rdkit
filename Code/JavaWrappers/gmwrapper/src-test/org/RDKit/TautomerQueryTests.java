
//
// Created by Gareth Jones on 6/1/2020.
//
// Copyright 2020 Schrodinger, Inc
//
package org.RDKit;

import static org.junit.Assert.*;
import org.junit.Test;

public class TautomerQueryTests extends  GraphMolTest {

    @Test
    public void basicOperations() {
        var mol = RWMol.MolFromSmiles("O=C1CCCCC1");
        var tautomerQuery = TautomerQuery.fromMol(mol);

        var target = RWMol.MolFromSmiles("OC1=CCCC(CC)C1");
        var match = tautomerQuery.isSubstructOf(target);
        assertTrue(match);
        var allMatches = tautomerQuery.substructOf(target);
        assertEquals(1, allMatches.size());
        var atomMatches = allMatches.get(0);
        assertEquals(7, atomMatches.size());

        var params = new SubstructMatchParameters();
        var matchingTautomers = new ROMol_Vect();
        allMatches = tautomerQuery.substructOf(target, params, matchingTautomers);
        assertEquals(1, allMatches.size());
        assertEquals(1, matchingTautomers.size());

        match = target.hasSubstructMatch(matchingTautomers.get(0));
        assertTrue(match);

        var tautomerFingerprint = tautomerQuery.patternFingerprintTemplate();
        var targetFingerprint = TautomerQuery.patternFingerprintTarget(target);
        var matchingFingerprints =  RDKFuncs.AllProbeBitsMatch(tautomerFingerprint, targetFingerprint);
        assertTrue(matchingFingerprints);

        tautomerFingerprint.delete();
        targetFingerprint.delete();
    }

}
