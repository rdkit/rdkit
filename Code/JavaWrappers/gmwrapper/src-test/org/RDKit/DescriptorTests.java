/* 
 * $Id: DescriptorTests.java 131 2011-01-20 22:01:29Z ebakke $
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

import org.junit.*;

public class DescriptorTests extends GraphMolTest {
    private File dataFile;

    public class DescriptorTestEntry {
	String smi;
	double clogp;
	double mr;
	int lineNo;

	public DescriptorTestEntry(int lineNo, String[] line) {
	    this.lineNo = lineNo;
	    smi = line[0].trim();
	    clogp = Double.parseDouble(line[1].trim());
	    mr = Double.parseDouble(line[2].trim());
	}
    }

    public List<DescriptorTestEntry> readData(String fn) throws Exception {
	List<DescriptorTestEntry> data = new ArrayList<DescriptorTestEntry>();
	int lineNo = 0;
	BufferedReader reader = new BufferedReader(new FileReader(fn));
	String line;
	// Skip the header
	try {
	    reader.readLine();
	    while ((line = reader.readLine()) != null) {
		lineNo += 1;
		if (line.length() > 0 && !line.startsWith("#")) {
		    String[] splitL = line.split("\t");
		    if (splitL.length == 3)
			data.add(new DescriptorTestEntry(lineNo, splitL));

		}
	    }
	} finally {
	    reader.close();
	}
	return data;
    }

    @Before
	public void setUp() {
	File base = getRdBase();
	dataFile = new File(base, "rdkit" + File.separator + "Chem" + File.separator + "test_data"
			    + File.separator + "Crippen.csv");
    }

    @Test
	public void testLogP() throws Exception {
	List<DescriptorTestEntry> data = readData(dataFile.getPath());
	for (DescriptorTestEntry item : data) {
	    ROMol mol = RWMol.MolFromSmiles(item.smi);
	    Double_Pair calcProps = RDKFuncs.calcCrippenDescriptors(mol);
	    assertEquals("cLogp on " + item.smi, item.clogp, calcProps.getFirst(), defaultDoubleTol);
	    assertEquals("mr on " + item.smi, item.mr, calcProps.getSecond(), defaultDoubleTol);

	}
    }

    @Test
	public void testRepeat() throws Exception {
	List<DescriptorTestEntry> data = readData(dataFile.getPath());
	for (DescriptorTestEntry item : data) {
	    ROMol mol = RWMol.MolFromSmiles(item.smi);
	    Double_Pair calcProps = RDKFuncs.calcCrippenDescriptors(mol);
	    assertEquals("cLogp on " + item.smi, item.clogp, calcProps.getFirst(), defaultDoubleTol);
	    assertEquals("mr on " + item.smi, item.mr, calcProps.getSecond(), defaultDoubleTol);

	}
    }

    @Test
	public void testIssue80() {
	ROMol m = RWMol.MolFromSmiles("CCOC");
	double ref = RDKFuncs.calcCrippenDescriptors(m).getFirst();
	double probe = RDKFuncs.calcCrippenDescriptors(m).getFirst();
	assertEquals(ref, probe, 0.0);
    }

    @Test
	public void testIssue1749494() {
	ROMol m1 = RWMol.MolFromSmiles("[*]CC");
	double v = RDKFuncs.calcCrippenDescriptors(m1).getFirst();
	assertEquals( 0.9739,v,defaultDoubleTol);
    }

    @Test public void testDescriptors1() {
	ROMol m1,m2,m3;
	m1 = RWMol.MolFromSmiles("C1=CC=CC=C1");
	m2 = RWMol.MolFromSmiles("C1=CC=CC=N1");
	m3 = RWMol.MolFromSmiles("C1=CC=CC=C1CC(=O)O");
        long tmp;
        tmp=RDKFuncs.calcNumHBA(m1);
        assertEquals(tmp,0);
        tmp=RDKFuncs.calcNumHBA(m2);
        assertEquals(tmp,1);
        tmp=RDKFuncs.calcNumHBA(m3);
        assertEquals(tmp,1);
        tmp=RDKFuncs.calcNumHBD(m1);
        assertEquals(tmp,0);
        tmp=RDKFuncs.calcNumHBD(m2);
        assertEquals(tmp,0);
        tmp=RDKFuncs.calcNumHBD(m3);
        assertEquals(tmp,1);
        tmp=RDKFuncs.calcLipinskiHBA(m1);
        assertEquals(tmp,0);
        tmp=RDKFuncs.calcLipinskiHBA(m2);
        assertEquals(tmp,1);
        tmp=RDKFuncs.calcLipinskiHBA(m3);
        assertEquals(tmp,2);
        tmp=RDKFuncs.calcLipinskiHBD(m1);
        assertEquals(tmp,0);
        tmp=RDKFuncs.calcLipinskiHBD(m2);
        assertEquals(tmp,0);
        tmp=RDKFuncs.calcLipinskiHBD(m3);
        assertEquals(tmp,1);
    }
    @Test public void testDescriptors2() {
	ROMol m1,m2;
	m1 = RWMol.MolFromSmiles("C1=CC=CC=C1");
	m2 = RWMol.MolFromSmiles("C1=CC=CC=N1");
        double tmp;
        tmp=RDKFuncs.calcMolLogP(m1);
        assertEquals(tmp,1.687,0.001);
        tmp=RDKFuncs.calcMolLogP(m2);
        assertEquals(tmp,1.082,0.001);
        tmp=RDKFuncs.calcMolMR(m1);
        assertEquals(tmp,26.442,0.001);
        tmp=RDKFuncs.calcMolMR(m2);
        assertEquals(tmp,24.237,0.001);

    }
    @Test public void testDescriptors3() {
	ROMol m1,m2,m3;
	m1 = RWMol.MolFromSmiles("C1=CC=CC=C1");
	m2 = RWMol.MolFromSmiles("C1=CC=CC=N1");
	m3 = RWMol.MolFromSmiles("C1=CC=CC=C1CC(=O)O");
        double tmp;
        tmp=RDKFuncs.calcTPSA(m1);
        assertEquals(tmp,0,0.01);
        tmp=RDKFuncs.calcTPSA(m2);
        assertEquals(tmp,12.89,0.01);
        tmp=RDKFuncs.calcTPSA(m3);
        assertEquals(tmp,37.30,0.01);
    }
    @Test public void testVSADescriptors() {
	ROMol m1,m2,m3;
	m1 = RWMol.MolFromSmiles("CO");
	m2 = RWMol.MolFromSmiles("CCO");
        Double_Vect v;

        v = RDKFuncs.calcSlogP_VSA(m1);
        assertEquals(v.get(0),0.0,0.001);
        assertEquals(v.get(1),12.216,0.001);
        v = RDKFuncs.calcSlogP_VSA(m2);
        assertEquals(v.get(0),0.0,0.001);
        assertEquals(v.get(1),11.713,0.001);
        assertEquals(v.get(4),6.924,0.001);

        v = RDKFuncs.calcSMR_VSA(m1);
        assertEquals(v.get(0),5.106,0.001);
        assertEquals(v.get(1),0.0,0.001);
        assertEquals(v.get(5),7.110,0.001);
        v = RDKFuncs.calcSMR_VSA(m2);
        assertEquals(v.get(0),5.106,0.001);
        assertEquals(v.get(4),6.924,0.001);
        assertEquals(v.get(5),6.607,0.001);

        v = RDKFuncs.calcPEOE_VSA(m1);
        assertEquals(v.get(0),5.106,0.001);
        assertEquals(v.get(7),7.110,0.001);
        assertEquals(v.get(5),0.000,0.001);
        v = RDKFuncs.calcPEOE_VSA(m2);
        assertEquals(v.get(0),5.106,0.001);
        assertEquals(v.get(6),6.924,0.001);
        assertEquals(v.get(7),6.607,0.001);
    }

    @Test public void testMW() {
	ROMol m1;
	m1 = RWMol.MolFromSmiles("C");
        assertEquals(16.043,RDKFuncs.calcAMW(m1),0.001);
        assertEquals(16.031,RDKFuncs.calcExactMW(m1),0.001);
    }
    @Test public void testMolFormula() {
	ROMol m1;
	m1 = RWMol.MolFromSmiles("C([2H])([3H])O");
        assertEquals("CH4O",RDKFuncs.calcMolFormula(m1));
        assertEquals("CH2DTO",RDKFuncs.calcMolFormula(m1,true));
        assertEquals("CH2[2H][3H]O",RDKFuncs.calcMolFormula(m1,true,false));
    }
    @Test public void testMolFormula2() {
	ROMol m1;
	m1 = RWMol.MolFromSmiles("C[13C]([2H])([3H])O");
        assertEquals("C2H6O",RDKFuncs.calcMolFormula(m1));
        assertEquals("C[13C]H4DTO",RDKFuncs.calcMolFormula(m1,true));
        assertEquals("C[13C]H4[2H][3H]O",RDKFuncs.calcMolFormula(m1,true,false));
    }

    @Test public void testConnectivityDescriptors1() {
	ROMol m1;
	m1 = RWMol.MolFromSmiles("c1ccccc1O");
        double tmp;
        tmp=RDKFuncs.calcChi0v(m1);
        assertEquals(tmp,3.834,0.01);
        tmp=RDKFuncs.calcChi1v(m1);
        assertEquals(tmp,2.134,0.01);
        tmp=RDKFuncs.calcChi2v(m1);
        assertEquals(tmp,1.336,0.01);
        tmp=RDKFuncs.calcChi3v(m1);
        assertEquals(tmp,0.756,0.01);
        tmp=RDKFuncs.calcChi4v(m1);
        assertEquals(tmp,0.428,0.01);
        tmp=RDKFuncs.calcChiNv(m1,4);
        assertEquals(tmp,0.428,0.01);

        tmp=RDKFuncs.calcChi0n(m1);
        assertEquals(tmp,3.834,0.01);
        tmp=RDKFuncs.calcChi1n(m1);
        assertEquals(tmp,2.134,0.01);
        tmp=RDKFuncs.calcChi2n(m1);
        assertEquals(tmp,1.336,0.01);
        tmp=RDKFuncs.calcChi3n(m1);
        assertEquals(tmp,0.756,0.01);
        tmp=RDKFuncs.calcChi4n(m1);
        assertEquals(tmp,0.428,0.01);
        tmp=RDKFuncs.calcChiNn(m1,4);
        assertEquals(tmp,0.428,0.01);
    }
    @Test public void testConnectivityDescriptors2() {
	ROMol m1;
        double tmp;

	m1 = RWMol.MolFromSmiles("C12CC2C3CC13");
        tmp=RDKFuncs.calcKappa1(m1);
        assertEquals(tmp,2.344,0.01);

	m1 = RWMol.MolFromSmiles("CC(C)C1CCC(C)CCC1");
        tmp=RDKFuncs.calcKappa2(m1);
        assertEquals(tmp,4.133,0.01);

	m1 = RWMol.MolFromSmiles("CC(C)C1CCC(C)CCC1");
        tmp=RDKFuncs.calcKappa3(m1);
        assertEquals(tmp,2.844,0.01);
    }



    public static void main(String args[]) {
	org.junit.runner.JUnitCore.main("org.RDKit.DescriptorTests");
    }

}
