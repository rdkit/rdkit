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
	    KeyFromPropHolder keys = new KeyFromPropHolder();
	    
            // fpholder - can we create it...
            PatternHolder ph = new PatternHolder();
	    TautomerPatternHolder tph = new TautomerPatternHolder();

            // Now lets make some molecules
            mol = RWMol.MolFromSmiles("c1ccccc1");
            String smiles = mol.MolToSmiles();
            mol.setProp("_Name", "foo");
	    
            mh.addMol(mol);
            assertEquals(1, mh.size());
            cmh.addMol(mol);
            assertEquals(1, cmh.size());
            csmh.addMol(mol);
            assertEquals(1, csmh.size());
            ctsmh.addMol(mol);
            assertEquals(1, ctsmh.size());
            keys.addMol(mol);
	    assertEquals(1, keys.size());
	    
            m = mh.getMol(0);
            assertEquals(smiles, m.MolToSmiles());
            m = cmh.getMol(0);
            assertEquals(smiles, m.MolToSmiles());
            m = csmh.getMol(0);
            assertEquals(smiles, m.MolToSmiles());
            m = ctsmh.getMol(0);
            assertEquals(smiles, m.MolToSmiles());
	    assertEquals("foo", keys.getKey(0));
	    
            mol = RWMol.MolFromSmiles("CCN");
	    mol.setProp("_Name", "bar");
            smiles = mol.MolToSmiles();

            mh.addMol(mol);
            assertEquals(2, mh.size());
            cmh.addMol(mol);
            assertEquals(2, mh.size());
	    keys.addMol(mol);
	    assertEquals(2, keys.size());
	    
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
	    assertEquals("bar", keys.getKey(1));
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
	    assertEquals(1, lib.countMatches(mol));
	    assertTrue(lib.hasMatch(mol));

	    
            if(lib.canSerialize()) {
                byte pickle[] = lib.Serialize();
                assertTrue(pickle != null);
                lib2 = SubstructLibrary.Deserialize(pickle);
                assertFalse(lib2 == null);
                assertEquals(lib.size(), 1);
                matches = lib.getMatches(mol);
                assertEquals(1, matches.size());
		assertEquals(1, lib.countMatches(mol));
		assertTrue(lib.hasMatch(mol));
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

  	@Test 
	public void test4Basics() {
            MolHolder mh = new MolHolder();
            TautomerPatternHolder pat = new TautomerPatternHolder();
            assertEquals(0, mh.size());
            CachedMolHolder cmh = new CachedMolHolder();
            assertEquals(0, cmh.size());
            CachedSmilesMolHolder csmh = new CachedSmilesMolHolder();
            assertEquals(0, csmh.size());
            CachedTrustedSmilesMolHolder ctsmh = new CachedTrustedSmilesMolHolder();
            assertEquals(0, ctsmh.size());
            mol = RWMol.MolFromSmiles("c1ccccc1");
	    TautomerQuery tautomerQuery = TautomerQuery.fromMol(mol);
	    
            // mol holder
            SubstructLibrary lib = new SubstructLibrary(mh, pat);
            lib.addMol(mol);
            
            UInt_Vect matches = lib.getMatches(mol);
            assertEquals(1, matches.size());

	    matches = lib.getMatches(tautomerQuery);
            assertEquals(1, matches.size());
	    assertEquals(1, lib.countMatches(mol));
	    assertEquals(1, lib.countMatches(tautomerQuery));
	    
            // cached mol holder
            try {
              lib = new SubstructLibrary(cmh, pat);
              lib.addMol(mol);
              assertTrue(false); // shouldn't get here can't share patternholder and addmol
            } catch(Exception e) {
              // ok to get here
            }
            cmh = new CachedMolHolder();
            pat = new TautomerPatternHolder();
            lib = new SubstructLibrary(cmh, pat);
            lib.addMol(mol);
            matches = lib.getMatches(mol);
            assertEquals(1, matches.size());
	    matches = lib.getMatches(tautomerQuery);
            assertEquals(1, matches.size());
	    assertEquals(1, lib.countMatches(mol));
	    assertEquals(1, lib.countMatches(tautomerQuery));
	    
            // cached smiles mol holder
            pat = new TautomerPatternHolder();
            lib = new SubstructLibrary(csmh, pat);
            lib.addMol(mol);

            matches = lib.getMatches(mol);
            assertEquals(1, matches.size());
	    matches = lib.getMatches(tautomerQuery);
            assertEquals(1, matches.size());
	    assertEquals(1, lib.countMatches(mol));
	    assertEquals(1, lib.countMatches(tautomerQuery));
	    
            // cached trusted smiles mol holder
            pat = new TautomerPatternHolder();
            lib = new SubstructLibrary(ctsmh, pat);
            lib.addMol(mol);

            matches = lib.getMatches(mol);
            assertEquals(1, matches.size());
	    matches = lib.getMatches(tautomerQuery);
            assertEquals(1, matches.size());
	    assertEquals(1, lib.countMatches(mol));
	    assertEquals(1, lib.countMatches(tautomerQuery));
        }
    

  	@Test 
	public void test5Basics() {
            MolHolder mh = new MolHolder();
            MolHolder mh2 = new MolHolder();
            PatternHolder pat = new PatternHolder();
            assertEquals(0, mh.size());
            CachedMolHolder cmh = new CachedMolHolder();
            assertEquals(0, cmh.size());
            mol = RWMol.MolFromSmiles("c1ccccc1");
	    mol.setProp("_Name", "foo");
	    KeyFromPropHolder keys = new KeyFromPropHolder();
	    KeyFromPropHolder keys2 = new KeyFromPropHolder();

            // mol holder
            SubstructLibrary lib = new SubstructLibrary(mh, pat, keys);
            lib.addMol(mol);
            
            SubstructLibrary lib2 = new SubstructLibrary(mh2, keys2);
	    lib2.addMol(mol);

            UInt_Vect matches = lib.getMatches(mol);
            UInt_Vect matches2 = lib2.getMatches(mol);
            assertEquals(1, matches.size());
	    assertEquals(matches.get(0), matches.get(0));

	    Str_Vect ids = lib.getKeys().getKeys(matches);
	    Str_Vect ids2 = lib2.getKeys().getKeys(matches2);
	    assertEquals(ids.get(0), "foo");
	    assertEquals(ids2.get(0), "foo");
        }
  
	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.SubstructLibraryTests");
	}

}

