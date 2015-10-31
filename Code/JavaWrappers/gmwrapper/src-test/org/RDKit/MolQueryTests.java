/* 
 * $Id: BasicMoleculeTests.java 131 2011-01-20 22:01:29Z ebakke $
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

public class MolQueryTests extends GraphMolTest {
    private ROMol qmol_sma,qmol_smi;
    @Before public void setUp() {
        String smiles="C1CCC1*";
        qmol_sma = RWMol.MolFromSmarts(smiles);
        qmol_smi = RWMol.MolFromSmiles(smiles);
    }
    @Test public void testBasicsSma() {
        ROMol m1,m2;
        m1 = RWMol.MolFromSmiles("C1CCC1C");
        m2 = RWMol.MolFromSmiles("C1CC(C)C1C");
        assertTrue(m1.hasSubstructMatch(qmol_sma));
        assertTrue(m2.hasSubstructMatch(qmol_sma));
        ROMol aqmol = RDKFuncs.adjustQueryProperties(qmol_sma);
        assertTrue(m1.hasSubstructMatch(aqmol));
        assertFalse(m2.hasSubstructMatch(aqmol));
    }
    @Test public void testBasicsSmi() {
        ROMol m1,m2;
        m1 = RWMol.MolFromSmiles("C1CCC1C");
        m2 = RWMol.MolFromSmiles("C1CC(C)C1C");
        assertFalse(m1.hasSubstructMatch(qmol_smi));
        assertFalse(m2.hasSubstructMatch(qmol_smi));
        ROMol aqmol = RDKFuncs.adjustQueryProperties(qmol_smi);
        assertTrue(m1.hasSubstructMatch(aqmol));
        assertFalse(m2.hasSubstructMatch(aqmol));
    }
    public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("org.RDKit.BasicMoleculeTests");
    }
}
