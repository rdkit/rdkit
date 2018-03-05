/* 
 * $Id: AtomPairsTests.java 4395 2011-01-28 16:54:48Z landrgr1 $
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

import org.junit.*;

public class AtomPairsTests extends GraphMolTest {

	@Test
	public void testAtomCodes_1() {
		ROMol mol = RWMol.MolFromSmiles("C=C");
		Atom atom0 = mol.getAtomWithIdx(0);
		Atom atom1 = mol.getAtomWithIdx(1);
		assertNotNull(atom0);
		assertNotNull(atom1);
		assertEquals(RDKFuncs.getAtomCode(atom0), RDKFuncs.getAtomCode(atom1));

		long tgt = 1 | (1 | 1 << RDKFuncs.getNumPiBits()) << RDKFuncs.getNumBranchBits();
		assertEquals(tgt, RDKFuncs.getAtomCode(atom0));
		tgt = 1 << RDKFuncs.getNumBranchBits() | 1 << (RDKFuncs.getNumBranchBits() + RDKFuncs.getNumPiBits());
		assertEquals(tgt, RDKFuncs.getAtomCode(atom0, 1));
	}

	@Test
	public void testAtomCodes_2() {
		ROMol mol = RWMol.MolFromSmiles("C#CO");
		long tgt = 1 | 2 << RDKFuncs.getNumBranchBits() | 1 << (RDKFuncs.getNumBranchBits() + RDKFuncs.getNumPiBits());
		assertEquals(tgt, RDKFuncs.getAtomCode(mol.getAtomWithIdx(0)));
		tgt = 2 | 2 << RDKFuncs.getNumBranchBits() | 1 << (RDKFuncs.getNumBranchBits() + RDKFuncs.getNumPiBits());
		assertEquals(tgt, RDKFuncs.getAtomCode(mol.getAtomWithIdx(1)));
		tgt = 1 | 0 << RDKFuncs.getNumBranchBits() | 3 << (RDKFuncs.getNumBranchBits() + RDKFuncs.getNumPiBits());
		assertEquals(tgt, RDKFuncs.getAtomCode(mol.getAtomWithIdx(2)));
	}

	@Test
	public void testAtomCodes_3() {
		ROMol mol = RWMol.MolFromSmiles("CC(O)C(O)(O)C");
		long tgt = 1 | 0 << RDKFuncs.getNumBranchBits() | 1 << (RDKFuncs.getNumBranchBits() + RDKFuncs.getNumPiBits());
		assertEquals(tgt, RDKFuncs.getAtomCode(mol.getAtomWithIdx(1), 2));
		tgt = 2 | 0 << RDKFuncs.getNumBranchBits() | 1 << (RDKFuncs.getNumBranchBits() + RDKFuncs.getNumPiBits());
		assertEquals(tgt, RDKFuncs.getAtomCode(mol.getAtomWithIdx(3), 2));

	}

	@Test
	public void testAtomCodes_4() {
		ROMol mol = RWMol.MolFromSmiles("C=CC(=O)O");
		long tgt = 0 | 0 << RDKFuncs.getNumBranchBits() | 3 << (RDKFuncs.getNumBranchBits() + RDKFuncs.getNumPiBits());
		assertEquals(tgt, RDKFuncs.getAtomCode(mol.getAtomWithIdx(4), 1));
		tgt = 3 | 1 << RDKFuncs.getNumBranchBits() | 1 << (RDKFuncs.getNumBranchBits() + RDKFuncs.getNumPiBits());
		assertEquals(tgt, RDKFuncs.getAtomCode(mol.getAtomWithIdx(2)));
	}

	@Test
	public void testAtomPairs() {
		ROMol mol = RWMol.MolFromSmiles("CCCCC");
		SparseIntVect32 fp;
		long tgt, c1, c2, c3;

		c1 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(0));
		c2 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(1));
		c3 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(2));
		tgt = 1 | (Math.min(c1, c2) | Math.max(c1, c2) << RDKFuncs.getCodeSize()) << RDKFuncs.getNumPathBits();
		assertEquals(tgt, RDKFuncs.getAtomPairCode(c1, c2, 1));
		assertEquals(tgt, RDKFuncs.getAtomPairCode(c2, c1, 1));
		tgt = 2 | (Math.min(c1, c3) | Math.max(c1, c3) << RDKFuncs.getCodeSize()) << RDKFuncs.getNumPathBits();
		assertEquals(tgt, RDKFuncs.getAtomPairCode(c1, c3, 2));
		assertEquals(tgt, RDKFuncs.getAtomPairCode(c3, c1, 2));

		mol = RWMol.MolFromSmiles("CCC");
		fp = RDKFuncs.getAtomPairFingerprint(mol);
		assertEquals(3, fp.getTotalVal());
		assertEquals(2, fp.getNonzero().size());

		c1 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(0));
		c2 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(1));
		c3 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(2));
		assertEquals(2, fp.getVal((int) RDKFuncs.getAtomPairCode(c1, c2, 1)));
		assertEquals(1, fp.getVal((int) RDKFuncs.getAtomPairCode(c1, c3, 2)));

		mol = RWMol.MolFromSmiles("CC=O.Cl");
		fp = RDKFuncs.getAtomPairFingerprint(mol);
		assertEquals(3, fp.getTotalVal());
		assertEquals(3, fp.getNonzero().size());

		c1 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(0));
		c2 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(1));
		c3 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(2));
		assertEquals(1, fp.getVal((int) RDKFuncs.getAtomPairCode(c1, c2, 1)));
		assertEquals(1, fp.getVal((int) RDKFuncs.getAtomPairCode(c1, c2, 1)));
		assertEquals(1, fp.getVal((int) RDKFuncs.getAtomPairCode(c2, c3, 1)));
	}

	@Test
	public void testAtomPairs2() {
		ROMol mol;
		SparseIntVect32 fp;

		mol = RWMol.MolFromSmiles("CCC");
		fp = RDKFuncs.getAtomPairFingerprint(mol, 1, 2);
		assertEquals(3, fp.getTotalVal());
		assertEquals(2, fp.getNonzero().size());

		fp = RDKFuncs.getAtomPairFingerprint(mol, 2, 2);
		assertEquals(1, fp.getTotalVal());
		assertEquals(1, fp.getNonzero().size());
	}

	@Test
	public void testHashedAtomPairs() {
		ROMol mol = RWMol.MolFromSmiles("c1ccccc1");
		SparseIntVect32 fp1, fp2;
		fp1 = RDKFuncs.getHashedAtomPairFingerprint(mol);
		fp2 = RDKFuncs.getHashedAtomPairFingerprint(mol);
		assertEquals(1.0, RDKFuncs.DiceSimilarity(fp1, fp2), 0.0);
		assertTrue( fp1.eq(fp2) );

		mol = RWMol.MolFromSmiles("c1ccccn1");
		fp2 = RDKFuncs.getHashedAtomPairFingerprint(mol);
		assertEquals(0.0, RDKFuncs.DiceSimilarity(fp1, fp2), 1.0);

		mol = RWMol.MolFromSmiles("c1ccccc1");
		fp1 = RDKFuncs.getHashedAtomPairFingerprint(mol, 2048);
		fp2 = RDKFuncs.getHashedAtomPairFingerprint(mol, 2048, 1, 3);
		assertEquals(1.0, RDKFuncs.DiceSimilarity(fp1, fp2), 0.0);
		assertTrue( fp1.eq(fp2) );

		fp2 = RDKFuncs.getHashedAtomPairFingerprint(mol, 2048, 1, 2);
		assertEquals(0.0, RDKFuncs.DiceSimilarity(fp1, fp2), 1.0);
	}

	@Test
	public void testTorsions() {
		ROMol mol = RWMol.MolFromSmiles("CCCC");
		SparseIntVect64 fp;
		double tgt;
		long c1, c2, c3, c4;
		UInt_Vect codes = new UInt_Vect();

		mol = RWMol.MolFromSmiles("CCCC");
		c1 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(0)) - 1;
		c2 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(1)) - 2;
		c3 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(2)) - 2;
		c4 = RDKFuncs.getAtomCode(mol.getAtomWithIdx(3)) - 1;
		tgt = c1 | (c2 | (c3 | c4 << RDKFuncs.getCodeSize()) << RDKFuncs.getCodeSize()) << RDKFuncs.getCodeSize();

		codes.clear();
		codes.add(c1);
		codes.add(c2);
		codes.add(c3);
		codes.add(c4);
		assertEquals(tgt, RDKFuncs.getTopologicalTorsionCode(codes).doubleValue(), 0.0);

		fp = RDKFuncs.getTopologicalTorsionFingerprint(mol);
		assertEquals(1, fp.getTotalVal());
		assertEquals(1, fp.getNonzero().size());

		mol = RWMol.MolFromSmiles("CCCCO.Cl");
		fp = RDKFuncs.getTopologicalTorsionFingerprint(mol);
		assertEquals(2, fp.getTotalVal());
		assertEquals(2, fp.getNonzero().size());

		fp = RDKFuncs.getTopologicalTorsionFingerprint(mol, 3);
		assertEquals(3, fp.getTotalVal());
		assertEquals(3, fp.getNonzero().size());
	}

	@Test
	public void testHashedTorsions() {
		ROMol mol = RWMol.MolFromSmiles("c1ccccc1");
		SparseIntVect64 fp1, fp2;
		fp1 = RDKFuncs.getHashedTopologicalTorsionFingerprint(mol);
		fp2 = RDKFuncs.getHashedTopologicalTorsionFingerprint(mol);
		assertEquals(1.0, RDKFuncs.DiceSimilarity(fp1, fp2), 0.0);
		assertTrue( fp1.eq(fp2) );

		mol = RWMol.MolFromSmiles("c1ccccn1");
		fp2 = RDKFuncs.getHashedTopologicalTorsionFingerprint(mol);
		assertEquals(0.0, RDKFuncs.DiceSimilarity(fp1, fp2), 1.0);

		mol = RWMol.MolFromSmiles("c1ccccc1");
		fp1 = RDKFuncs.getHashedTopologicalTorsionFingerprint(mol, 2048, 6);
		fp2 = RDKFuncs.getHashedTopologicalTorsionFingerprint(mol, 2048, 6);
		assertEquals(1.0, RDKFuncs.DiceSimilarity(fp1, fp2), 0.0);
		assertTrue( fp1.eq(fp2) );

		mol = RWMol.MolFromSmiles("c1ccccn1");
		fp2 = RDKFuncs.getHashedTopologicalTorsionFingerprint(mol, 2048, 6);
		assertEquals(0.0, RDKFuncs.DiceSimilarity(fp1, fp2), 1.0);
	}

	@Test
	public void testBulkTorsions() {
		String fName = new File(getRdBase(), "/Projects/DbCLI/testData/pubchem.200.sdf").getPath();
		SDMolSupplier suppl = new SDMolSupplier(fName);
		while (!suppl.atEnd()) {
			ROMol mol = suppl.next();
			SparseIntVect64 fp;
			fp = RDKFuncs.getTopologicalTorsionFingerprint(mol);
			assertTrue(fp.getTotalVal() > 1);
		}
	}

	@Test
	public void testRootedAtomPairs() {
		ROMol mol = RWMol.MolFromSmiles("OCCCCC");
		SparseIntVect32 fp1, fp2;
		UInt_Vect roots = new UInt_Vect();

		fp1 = RDKFuncs.getAtomPairFingerprint(mol);
		Match_Vect nz1 = fp1.getNonzero();
		assertTrue(nz1.size() > 0);

		roots.add(0);
		fp2 = RDKFuncs.getAtomPairFingerprint(mol, roots);
		Match_Vect nz2 = fp2.getNonzero();
		assertTrue(nz2.size() > 0);
		assertTrue(nz2.size() < nz1.size());

		for (int i = 0; i < nz2.size(); i++) {
			Int_Pair pair = nz2.get(i);
			assertEquals(pair.getSecond() , fp2.getVal(pair.getFirst()));
		}
	}

	@Test
	public void testRootedTorsions() {
		SparseIntVect64 fp1, fp2;
		
		ROMol mol = RWMol.MolFromSmiles("OCCCC");
		UInt_Vect roots = new UInt_Vect();
		roots.add(0);
		
		fp1 = RDKFuncs.getTopologicalTorsionFingerprint(mol);
		Long_Pair_Vect nz1 = fp1.getNonzero();
		assertTrue(nz1.size() > 0);

		fp2 = RDKFuncs.getTopologicalTorsionFingerprint(mol, 4, roots);
		Long_Pair_Vect nz2 = fp2.getNonzero();
		assertTrue(nz2.size() > 0);
		assertTrue(nz2.size() < nz1.size());

		for (int i = 0; i < nz2.size(); ++i) {
			Long_Pair pair = nz2.get(i);
			assertTrue(pair.getSecond() <= fp2.getVal(pair.getFirst()));
		}
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main(AtomPairsTests.class.getName());
	}
}
