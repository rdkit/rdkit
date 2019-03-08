/* 
 *
 *  Copyright (c) 2018, Novartis Institutes for BioMedical Research Inc.
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
    
public class SubstructLibraryTests extends GraphMolTest {
	private ROMol mol;
        private ROMol m;

	@Before public void setUp() {
	}

	@Test 
	public void test1Basics() {
            SubstructLibrary lib = new SubstructLibrary();
            assertEquals(0, lib.size());
            // mol holders
            MolHolder mh = new MolHolder();
            assertEquals(0, mh.size());
            CachedMolHolder cmh = new CachedMolHolder();
            assertEquals(0, cmh.size());
            CachedSmilesMolHolder csmh = new CachedSmilesMolHolder();
            assertEquals(0, csmh.size());
            CachedTrustedSmilesMolHolder ctsmh = new CachedTrustedSmilesMolHolder();
            assertEquals(0, ctsmh.size());

            // fpholder - can we create it...
            PatternHolder ph = new PatternHolder();

            // Now lets make some molecules
            mol = RWMol.MolFromSmiles("c1ccccc1");
            String smiles = mol.MolToSmiles();
            
            mh.addMol(mol);
            assertEquals(1, mh.size());
            cmh.addMol(mol);
            assertEquals(1, cmh.size());
            csmh.addMol(mol);
            assertEquals(1, csmh.size());
            ctsmh.addMol(mol);
            assertEquals(1, ctsmh.size());
            
            m = mh.getMol(0);
            assertEquals(smiles, m.MolToSmiles());
            m = cmh.getMol(0);
            assertEquals(smiles, m.MolToSmiles());
            m = csmh.getMol(0);
            assertEquals(smiles, m.MolToSmiles());
            m = ctsmh.getMol(0);
            assertEquals(smiles, m.MolToSmiles());

            mol = RWMol.MolFromSmiles("CCN");
            smiles = mol.MolToSmiles();

            mh.addMol(mol);
            assertEquals(2, mh.size());
            cmh.addMol(mol);
            assertEquals(2, mh.size());
            
            csmh.addSmiles("CCN");
            assertEquals(2, csmh.size());
            ctsmh.addSmiles("CCN");
            assertEquals(2, ctsmh.size());

            m = mh.getMol(1);
            assertEquals(smiles, m.MolToSmiles());
            m = cmh.getMol(1);
            assertEquals(smiles, m.MolToSmiles());
            m = csmh.getMol(1);
            assertEquals(smiles, m.MolToSmiles());
            m = ctsmh.getMol(1);
            assertEquals(smiles, m.MolToSmiles());
	}


  	@Test 
	public void test2Basics() {
            MolHolder mh = new MolHolder();
            assertEquals(0, mh.size());
            CachedMolHolder cmh = new CachedMolHolder();
            assertEquals(0, cmh.size());
            CachedSmilesMolHolder csmh = new CachedSmilesMolHolder();
            assertEquals(0, csmh.size());
            CachedTrustedSmilesMolHolder ctsmh = new CachedTrustedSmilesMolHolder();
            assertEquals(0, ctsmh.size());
            mol = RWMol.MolFromSmiles("c1ccccc1");
            // mol holder
            SubstructLibrary lib = new SubstructLibrary(mh);
            SubstructLibrary lib2;
            lib.addMol(mol);
            
            UInt_Vect matches = lib.getMatches(mol);
            assertEquals(1, matches.size());

            if(lib.canSerialize()) {
                byte pickle[] = lib.Serialize();
                assertTrue(pickle != null);
                lib2 = SubstructLibrary.Deserialize(pickle);
                assertFalse(lib2 == null);
                assertEquals(lib.size(), 1);
                matches = lib.getMatches(mol);
                assertEquals(1, matches.size());
            }
            
            // cached mol holder
            lib = new SubstructLibrary(cmh);
            lib.addMol(mol);

            matches = lib.getMatches(mol);
            assertEquals(1, matches.size());

            if(lib.canSerialize()) {
                byte pickle[] = lib.Serialize();
                assertTrue(pickle != null);
                lib2 = SubstructLibrary.Deserialize(pickle);
                assertFalse(lib2 == null);
                assertEquals(lib.size(), 1);
                matches = lib.getMatches(mol);
                assertEquals(1, matches.size());
            }

            // cached smiles mol holder
            lib = new SubstructLibrary(csmh);
            lib.addMol(mol);

            matches = lib.getMatches(mol);
            assertEquals(1, matches.size());

            if(lib.canSerialize()) {
                byte pickle[] = lib.Serialize();
                assertTrue(pickle != null);
                lib2 = SubstructLibrary.Deserialize(pickle);
                assertFalse(lib2 == null);
                assertEquals(lib.size(), 1);
                matches = lib.getMatches(mol);
                assertEquals(1, matches.size());
            }
            
            // cached trusted smiles mol holder
            lib = new SubstructLibrary(ctsmh);
            lib.addMol(mol);

            matches = lib.getMatches(mol);
            assertEquals(1, matches.size());

            if(lib.canSerialize()) {
                byte pickle[] = lib.Serialize();
                assertTrue(pickle != null);
                lib2 = SubstructLibrary.Deserialize(pickle);
                assertFalse(lib2 == null);
                assertEquals(lib.size(), 1);
                matches = lib.getMatches(mol);
                assertEquals(1, matches.size());
            }
            
        }

  	@Test 
	public void test3Basics() {
            MolHolder mh = new MolHolder();
            PatternHolder pat = new PatternHolder();
            assertEquals(0, mh.size());
            CachedMolHolder cmh = new CachedMolHolder();
            assertEquals(0, cmh.size());
            CachedSmilesMolHolder csmh = new CachedSmilesMolHolder();
            assertEquals(0, csmh.size());
            CachedTrustedSmilesMolHolder ctsmh = new CachedTrustedSmilesMolHolder();
            assertEquals(0, ctsmh.size());
            mol = RWMol.MolFromSmiles("c1ccccc1");
            // mol holder
            SubstructLibrary lib = new SubstructLibrary(mh, pat);
            lib.addMol(mol);
            
            UInt_Vect matches = lib.getMatches(mol);
            assertEquals(1, matches.size());

            // cached mol holder
            try {
              lib = new SubstructLibrary(cmh, pat);
              lib.addMol(mol);
              assertTrue(false); // shouldn't get here can't share patternholder and addmol
            } catch(Exception e) {
              // ok to get here
            }
            cmh = new CachedMolHolder();
            pat = new PatternHolder();
            lib = new SubstructLibrary(cmh, pat);
            lib.addMol(mol);
            matches = lib.getMatches(mol);
            assertEquals(1, matches.size());

            // cached smiles mol holder
            pat = new PatternHolder();
            lib = new SubstructLibrary(csmh, pat);
            lib.addMol(mol);

            matches = lib.getMatches(mol);
            assertEquals(1, matches.size());
            
            // cached trusted smiles mol holder
            pat = new PatternHolder();
            lib = new SubstructLibrary(ctsmh, pat);
            lib.addMol(mol);

            matches = lib.getMatches(mol);
            assertEquals(1, matches.size());
        }

    
  
	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.SubstructLibraryTests");
	}

}
