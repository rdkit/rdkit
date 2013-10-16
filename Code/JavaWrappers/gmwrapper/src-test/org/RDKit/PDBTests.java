/* 
 * $Id: ForceFieldsTests.java 131 2011-01-20 22:01:29Z ebakke $
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

import java.io.File;

import org.junit.Test;

public class PDBTests extends GraphMolTest {

	File testDataDir = new File(getRdBase(), 
	"Code/GraphMol/FileParsers/test_data");

	@Test
	public void testPDBReadBasic () {
            File molFile = new File(testDataDir, "1CRN.pdb");
            ROMol m = RWMol.MolFromPDBFile(molFile.getPath());
            assertEquals(327, m.getNumAtoms());
            assertEquals(337, m.getNumBonds());
            AtomMonomerInfo mi=new AtomMonomerInfo(m.getAtomWithIdx(0).getMonomerInfo());
            assert(mi instanceof AtomPDBResidueInfo);
            // FIX: need to actually test this
	}

	@Test
	public void testPDBWriteBasic () {
            File molFile = new File(testDataDir, "1CRN.pdb");
            ROMol m = RWMol.MolFromPDBFile(molFile.getPath());
            String mb = new String(m.MolToPDBBlock());
            ROMol m2 = RWMol.MolFromPDBBlock(mb);
            assertEquals(327, m2.getNumAtoms());
            assertEquals(337, m2.getNumBonds());
            AtomMonomerInfo mi=new AtomMonomerInfo(m2.getAtomWithIdx(0).getMonomerInfo());
            assert(mi instanceof AtomPDBResidueInfo);
            // FIX: need to actually test this
            
	}


	public static void main(String args[]) {
            org.junit.runner.JUnitCore.main("org.RDKit.PDBTests");
	}

}
