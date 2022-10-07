/* 
 * $Id: WrapperTests.java 131 2011-01-20 22:01:29Z ebakke $
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
import java.io.*;
import java.io.File;
import java.util.ArrayList;

import org.junit.*;


public class WrapperTests extends GraphMolTest {

	private ArrayList<String> tmpFiles = new ArrayList<String>();

	@After
	public void tearDown() {
		for (String fn : tmpFiles) {
			new File(fn).delete();
		}
	}
	@Test 
    public void testBasicInstantiation_Conformer() {
        Conformer conformer = new Conformer();
        assertNotNull(conformer);
    }
    @Test 
    public void testBasicInstantiation_QueryAtom() {
	QueryAtom queryAtom = new QueryAtom();
        assertNotNull(queryAtom);
    }
    @Test
    public void testBasicInstantiation_QueryBond() {
	QueryBond queryBond = new QueryBond();
        assertNotNull(queryBond);
    }
    @Test
    public void testBasicInstantiation_QueryOps() {
        // Comes from QueryOps.cpp
        RecursiveStructureQuery rsQuery = new RecursiveStructureQuery();
        assertNotNull(rsQuery);

        AtomRingQuery arQuery = new AtomRingQuery();
        assertNotNull(arQuery);
    }
    @Test 
    public void testBasicInstantiation_Atom() {
	Atom atom = new Atom();
	assertNotNull(atom);
    }
    @Test 
    public void testBasicInstantiation_ROMol() {
	ROMol roMol = new ROMol();
	assertNotNull(roMol);
    }
    @Test 
    public void testBasicInstantiation_Bond() {
	Bond bond = new Bond();
	assertNotNull(bond);
    }
    @Test 
    public void testBasicInstantiation_MolStackElem() {
	// Came from Canon.cpp
	MolStackElem molStack = new MolStackElem(1);
	assertNotNull(molStack);
    }
    @Test 
    public void testBasicInstantiation_RWMol() {
	RWMol rwMol = new RWMol();
	assertNotNull(rwMol);
    }
    @Test
    public void testBasicInstantiation_PeriodicTable() {
        assertNotNull(PeriodicTable.getTable());
    }
    @Test
    public void testBasicInstantiation_MolSanitizeException() {
        MolSanitizeException sanit = new MolSanitizeException("some error message");
        assertNotNull(sanit);
    }
    @Test
    public void testBasicInstantiation_SmilesParseException() {
        SmilesParseException parse = new SmilesParseException("some error message");
        assertNotNull(parse);
    }
    @Test
    public void testBasicInstantiation_RingInfo() {
        RingInfo rInfo = new RingInfo();
        assertNotNull(rInfo);
    }
    @Test
    public void testBasicInstantiation_ChemicalReaction() {
        ChemicalReaction reaction = new ChemicalReaction();
        assertNotNull(reaction);
    }
       
    @Test
    public void testBasicInstantiation_BondIterator() {
       	ROMol m = RWMol.MolFromSmiles("CS");
    	BondIterator bonditer = new BondIterator(m);
        assertNotNull(bonditer);
    }
    @Test
    public void testBasicInstantiation_ConstBondIterator() {
       	ROMol m = RWMol.MolFromSmiles("CS");
    	ConstBondIterator bonditer = new ConstBondIterator(m);
        assertNotNull(bonditer);
    }
    @Test
    public void testBasicInstantiation_AtomIterator() {
       	ROMol m = RWMol.MolFromSmiles("CS");
        AtomIterator atomiter = new AtomIterator(m);
        assertNotNull(atomiter);
    }
    @Test
    public void testBasicInstantiation_HeteroatomIterator() {
    	ROMol m = RWMol.MolFromSmiles("CS");
    	HeteroatomIterator atomiter = new HeteroatomIterator(m);
        assertNotNull(atomiter);
    }
    @Test
    public void testBasicInstantiation_QueryAtomIterator() {
       	ROMol m = RWMol.MolFromSmiles("CS");
    	QueryAtomIterator atomiter = new QueryAtomIterator(m, new QueryAtom(6));
        assertNotNull(atomiter);
    }
    @Test
    public void testBasicInstantiation_AromaticAtomIterator() {
       	ROMol m = RWMol.MolFromSmiles("Cc1ccccc1");
    	AromaticAtomIterator atomiter = new AromaticAtomIterator(m);
        assertNotNull(atomiter);
    }
    @Test
    public void testBasicInstantiation_SDMolSupplier() {
        SDMolSupplier sup = new SDMolSupplier();
        assertNotNull(sup);
    }
    @Test
    public void testBasicInstantiation_ForwardSDMolSupplier() {
    	ForwardSDMolSupplier sup = new ForwardSDMolSupplier();
        assertNotNull(sup);
    }
    @Test
    public void testBasicInstantiation_SmilesMolSupplier() {
        SmilesMolSupplier sup = new SmilesMolSupplier();
        assertNotNull(sup);
    }
    @Test
    public void testBasicInstantiation_TDTMolSupplier() {
        TDTMolSupplier sup = new TDTMolSupplier();
        assertNotNull(sup);
    }

    @Test
    public void testBasicInstantiation_ResonanceMolSupplier() {
       	ROMol m = RWMol.MolFromSmiles("CC(=O)[O-]");
        ResonanceMolSupplier sup = new ResonanceMolSupplier(m);
        assertNotNull(sup);
    }

    @Test
    public void testBasicInstantiation_SDMWriter() {
        SDWriter mw = new SDWriter("tmp.sdf");
        mw.close();
        tmpFiles.add("tmp.sdf");
        assertNotNull(mw);
    }

    @Test
    public void testBasicInstantiation_TDTWriter() {
        TDTWriter mw = new TDTWriter("tmp.tdt");
        mw.close();
        tmpFiles.add("tmp.tdt");
        assertNotNull(mw);
    }
    @Test
    public void testBasicInstantiation_SmilesWriter() {
        SmilesWriter mw = new SmilesWriter("tmp.smi");
        mw.close();
        tmpFiles.add("tmp.smi");
        assertNotNull(mw);
    }
    @Test
    public void testBasicInstantiation_Transform2D() {
        Transform2D tr = new Transform2D();
        assertNotNull(tr);
    }
    @Test
    public void testBasicInstantiation_Transform3D() {
        Transform3D tr = new Transform3D();
        assertNotNull(tr);
    }

    //Tests for Exception#getMessage() wrapping
    private static final String TEST_MESSAGE = "test message";

    @Test
    public void testChemicalReactionException() {
        final ChemicalReactionException e =
                new ChemicalReactionException(TEST_MESSAGE);
        assertEquals(TEST_MESSAGE, e.getMessage());

    }

    @Test
    public void testChemicalReactionParserException() {
        final ChemicalReactionParserException e =
                new ChemicalReactionParserException(TEST_MESSAGE);
        assertEquals(TEST_MESSAGE, e.getMessage());

    }

    @Test
    public void testConformerException() {
        ROMol mol = null;
        try {
            mol = RWMol.MolFromSmiles("c1ccccc1");
            mol.getConformer();
        } catch (ConformerException e) {
            assertEquals("No conformations available on the molecule",
                    e.getMessage());
        } finally {
            if (mol != null) {
                mol.delete();
            }
        }

    }

    @Test
    public void testMolPicklerException() {
        final MolPicklerException e = new MolPicklerException(TEST_MESSAGE);
        assertEquals(TEST_MESSAGE, e.getMessage());

    }

    @Test
    public void testMolSanitizeException() {
        RWMol mol = null;
        try {
            mol = RWMol.MolFromSmiles("c1cccc1");
            mol.sanitizeMol();
        } catch (MolSanitizeException e) {
            assertEquals("Can't kekulize mol.  Unkekulized atoms: 0 1 2 3 4",
                    e.getMessage());
        } finally {
            if (mol != null) {
                mol.delete();
            }
        }
    }

    @Test
    public void testSmilesParseException() {
        SmilesParseException e = new SmilesParseException(TEST_MESSAGE);
        assertEquals(TEST_MESSAGE, e.getMessage());
    }

    @Test
    public void testGenericRDKitException() {
        final GenericRDKitException e =
                new GenericRDKitException(TEST_MESSAGE);
        assertEquals(TEST_MESSAGE, e.getMessage());

    }

    @Test
    public void testAtomSanitizeException() {
        final AtomSanitizeException e =
                new AtomSanitizeException(TEST_MESSAGE, 1L);
        assertEquals(TEST_MESSAGE, e.getMessage());

    }

    @Test
    public void testAtomValenceException() {
        final AtomValenceException e =
                new AtomValenceException(TEST_MESSAGE, 1L);
        assertEquals(TEST_MESSAGE, e.getMessage());
    }

    @Test
    public void testAtomKekulizeException() {
        final AtomKekulizeException e =
                new AtomKekulizeException(TEST_MESSAGE, 1L);
        assertEquals(TEST_MESSAGE, e.getMessage());
    }

    @Test
    public void testKekulizeException() {
        final UInt_Vect uInt_Vect = new UInt_Vect();
        uInt_Vect.add(1L);
        final KekulizeException e =
                new KekulizeException(TEST_MESSAGE, uInt_Vect);
        assertEquals(TEST_MESSAGE, e.getMessage());
        uInt_Vect.delete();
        e.delete();
    }

    @Test
    public void cdxmlReader() {
	String rdpath = System.getenv("RDBASE");
	if (rdpath == null)
	    org.junit.Assert.fail("No definition for RDBASE");
	File base = new File(rdpath);
	File testFile = new File(base, "Code" + File.separator + "GraphMol"
				 + File.separator + "test_data" + File.separator +
				 "CDXML" + File.separator + "beta-cypermethrin.cdxml");
	String fn = testFile.getAbsolutePath();
	RWMol_Vect prods = RWMol.MolsFromCDXMLFile(fn);
	assertEquals(prods.size(), 1);
	for(int idx = 0; idx < prods.size(); idx++) {
	    if(idx == 0) {
		assertEquals(prods.get(idx).MolToSmiles(true), "CC1(C)[C@H](C=C(Cl)Cl)[C@H]1C(=O)O[C@@H](C#N)c1cccc(Oc2ccccc2)c1");
	    }
	}

    }
    
    public static void main(String args[]) {
      org.junit.runner.JUnitCore.main("org.RDKit.WrapperTests");
    }
}


