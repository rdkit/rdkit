/* 
 * $Id: AromaticTests.java 131 2011-01-20 22:01:29Z ebakke $
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
import java.util.*;

import org.junit.Test;

public class AromaticTests extends GraphMolTest {

	// Check all lines in all files
	@Test
	public void testAromaticRegression() throws Exception {
		String filePath = getFilePath("aromat_regress.txt");
		performAromaticTest(filePath,0);
	}

	@Test
	public void testNCIAromaticRegression() throws Exception {
		String filePath = getFilePath("NCI_aromat_regress.txt");
		// assertEquals(filePath,"foo");
		performAromaticTest(filePath,0);
	}

	public void performAromaticTest(String filePath, int expectedFailures) throws Exception {
		List<AromaticTestEntry> testData = readData(filePath);
		int nFailed = 0;
		for (AromaticTestEntry test : testData) {
			ROMol mol = RWMol.MolFromSmiles(test.smi1);
			assertNotNull(mol);
			if (!test.smi2.equals("FAIL")) {
				int count = 0;
				ArrayList<Atom> aroms = new ArrayList<Atom>();
				Atom_Vect atoms = mol.getAtoms();
				for (int i = 0; i < atoms.size(); i++)
					if (atoms.get(i).getIsAromatic()) {
						count++;
						aroms.add(atoms.get(i));
					}
				if (test.numAromatics != count) {
					nFailed++;
				}
			}
		}
		assertEquals("More than " + expectedFailures + " on file " + filePath, expectedFailures, nFailed);

	}

	public String getFilePath(String fileName) {
		File base = getRdBase();
		File testFileDir = new File(base, "rdkit" + File.separator + "Chem" + File.separator + "test_data");
		return testFileDir.getAbsolutePath() + File.separator + fileName;
	}

	public class AromaticTestEntry {
		String smi1;
		String smi2;
		int numAromatics;
		int[] aromaticIdx; 
		int lineNo;

		public AromaticTestEntry(int lineNo, String[] line) {
			this.lineNo = lineNo;
			smi1 = line[0];
			smi2 = line[1];
			numAromatics = Integer.parseInt(line[2]);
			String aromaticList = line[3].replace("[", "").replace("]", "").trim();
			// Catch a null string
			String[] aromatics = aromaticList.length() > 0 ? aromaticList.split(",")
					: new String[0];
			assertTrue("bad test line at " + lineNo, aromatics.length == numAromatics || smi2.equals("FAIL"));
			aromaticIdx = new int[aromatics.length];
			for (int i = 0; i < aromatics.length; i++)
				aromaticIdx[i] = Integer.parseInt(aromatics[i].trim());
		}
	}

	public List<AromaticTestEntry> readData(String fn) throws Exception {
		List<AromaticTestEntry> data = new ArrayList<AromaticTestEntry>();
		int lineNo = 0;

		BufferedReader reader = new BufferedReader(new FileReader(fn));
		String line;
		// Skip the header
		try {
			reader.readLine();
			while ((line = reader.readLine()) != null) {
				lineNo++;
				if (line.length() > 0 && !line.startsWith("#")) {
					String[] splitL = line.split("\t");
					if (splitL.length == 4)
						data.add(new AromaticTestEntry(lineNo, splitL));

				}
			}
		} finally {
			reader.close();
		}
		return data;
	}

	public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.AromaticTests");
	}

}
