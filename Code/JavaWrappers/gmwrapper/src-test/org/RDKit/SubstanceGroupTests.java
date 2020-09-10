//
// Created by Gareth Jones on 9/4/2020.
//
// Copyright 2020 Schrodinger, Inc
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

package org.RDKit;

import org.junit.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

public class SubstanceGroupTests extends GraphMolTest {

    private RWMol getMol(String file) {
        file = getRdBase() + "/Code/GraphMol/test_data/"+file;
        try {
            String block = new String(Files.readAllBytes(Paths.get(file)));
            return RDKFuncs.MolBlockToMol(block);
        } catch (IOException ex) {
            throw new RuntimeException("IO Exception reading " + file, ex);
        }
    }

    @Test
    public void supGroupTest() {
        RWMol mol = getMol("sgroups_and_marvin_metal.mol");
        assertNotNull(mol);
        long sgroupCount = RDKFuncs.getSubstanceGroupCount(mol);
        assertEquals(1, sgroupCount);
        SubstanceGroup supGroup = RDKFuncs.getSubstanceGroupWithIdx(mol, 0);
        assertTrue(supGroup.hasProp("TYPE"));
        assertTrue(supGroup.hasProp("index"));
        assertTrue(supGroup.hasProp("LABEL"));
        assertTrue(supGroup.hasProp("DATAFIELDS"));
        assertEquals("SUP", supGroup.getStringProp("TYPE"));
        assertEquals(1, supGroup.getUIntProp("index"));
        assertEquals("M", supGroup.getStringProp("LABEL"));
        assertEquals(0, supGroup.getStringVectProp("DATAFIELDS").size());
        assertEquals(1, supGroup.getAtoms().size());
        assertEquals(0, supGroup.getBonds().size());
        assertEquals(1, supGroup.getAtoms().get(0));
    }

    @Test
    public void polymerTest() {

        RWMol mol = getMol("sgroups_copolymer2.mol");
        assertNotNull(mol);
        long sgroupCount = RDKFuncs.getSubstanceGroupCount(mol);
        assertEquals(3, sgroupCount);
        assertEquals(9, mol.getNumAtoms());

        SubstanceGroup substanceGroup = RDKFuncs.getSubstanceGroupWithIdx(mol, 0);
        assertTrue(substanceGroup.hasProp("index"));
        assertEquals(2, substanceGroup.getUIntProp("index") );
        assertEquals(2, substanceGroup.getAtoms().size());
        assertEquals(3, substanceGroup.getAtoms().get(0));
        assertEquals(2, substanceGroup.getAtoms().get(1));
        assertEquals(2, substanceGroup.getBonds().size());;
        assertEquals(1, substanceGroup.getBonds().get(0));
        assertEquals(3, substanceGroup.getBonds().get(1));
        assertTrue(substanceGroup.hasProp("PARENT"));
        assertEquals(10, substanceGroup.getUIntProp("PARENT"));

        substanceGroup = RDKFuncs.getSubstanceGroupWithIdx(mol, 2);
        assertTrue(substanceGroup.hasProp("index"));
        assertEquals(10, substanceGroup.getUIntProp("index") );
        assertEquals(5, substanceGroup.getAtoms().size());
        assertEquals(3, substanceGroup.getAtoms().get(0));
        assertEquals(2, substanceGroup.getAtoms().get(1));
        assertEquals(4, substanceGroup.getAtoms().get(2));
        assertEquals(5, substanceGroup.getAtoms().get(3));
        assertEquals(7, substanceGroup.getAtoms().get(4));
        assertEquals(2, substanceGroup.getBonds().size());;
        assertEquals(1, substanceGroup.getBonds().get(0));
        assertEquals(5, substanceGroup.getBonds().get(1));
        assertFalse(substanceGroup.hasProp("PARENT"));

    }


    public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("org.RDKit.SubstanceGroupTests");
    }

}
