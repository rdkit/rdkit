/*
* $Id: Chemv2Tests.java 131 2011-01-20 22:01:29Z ebakke $
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
*       products derived from this software without specific prior written
* permission.
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
import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedReader;
import java.util.Arrays;

import org.junit.Test;

public class Chemv2Tests extends GraphMolTest {

    /* Pickling tests skipped for the time being */
    @Test
    public void testBasicStuff() {
        ROMol m = RWMol.MolFromSmiles("COC(=O)O");
        Atom a1 = m.getAtomWithIdx(1);
        assertEquals( 8,a1.getAtomicNum() );
        assertEquals( 6,m.getAtomWithIdx(2).getAtomicNum() );
        Bond b1 = m.getBondWithIdx(1);
        assertEquals( Bond.BondType.SINGLE,b1.getBondType() );
        assertEquals( Bond.BondType.DOUBLE,m.getBondWithIdx(2).getBondType() );
        assertEquals( Bond.BondType.SINGLE,m.getBondBetweenAtoms(0, 1).getBondType() );
    }

    @Test
    public void testEditingPersisting()    {
        RWMol m = RWMol.MolFromSmiles("COC(=C)O");
        Atom a1 = m.getAtomWithIdx(3);
        assertEquals("bad atom order",6, a1.getAtomicNum());
        a1.setAtomicNum(7);
        assertEquals("bad atom order",7, a1.getAtomicNum());
        assertEquals("atom order not stored",7, m.getAtomWithIdx(3).getAtomicNum());
    }

    @Test
    public void testSMARTSBasics () {
        ROMol m = RWMol.MolFromSmiles("COC(=O)O");
        ROMol p = RWMol.MolFromSmarts("CO");
        assertTrue(m.hasSubstructMatch(p));
        ROMol p2 = RWMol.MolFromSmarts("CS");
        assertFalse(m.hasSubstructMatch(p2));
        assertEquals( 2,p.getNumAtoms() );
        assertEquals( 1,p.getNumBonds() );
        assertTrue(m.hasSubstructMatch(p));
        Match_Vect_Vect matches = m.getSubstructMatches(p);
        assertEquals( 3,matches.size() );
        Match_Vect match = matches.get(0);
        assertEquals("bad match length", 2, match.size() );
        matches = m.getSubstructMatches(p, false);
        assertEquals( 3,matches.size() );
        match = matches.get(0);
        assertEquals("bad match length", 2, match.size() );

        p = RWMol.MolFromSmarts("COC");
        assertTrue(m.hasSubstructMatch(p));
        assertEquals( 3,p.getNumAtoms() );
        assertEquals( 2,p.getNumBonds() );
        assertTrue(m.hasSubstructMatch(p));
        matches = m.getSubstructMatches(p);
        assertEquals( 1,matches.size() );
        matches = m.getSubstructMatches(p, false);
        assertEquals( 2,matches.size() );
    }

    @Test
    public void testDataGetSetSuccess() {
        ROMol m = RWMol.MolFromSmiles("CCOC");
        m.setProp("foo", "3");
        String v = m.getProp("foo");
        assertEquals("3",v);
    }

        @Test(expected=KeyErrorException.class)
    public void testDataGetSetFailure() {
        ROMol m = RWMol.MolFromSmiles("CCOC");
        m.getProp("monkey");
    }

    @Test
    public void testIssue399() {
        ROMol m = RWMol.MolFromSmiles("[C@H]1(C)CO1");
        m.compute2DCoords();
        Conformer c = m.getConformer();
        m.WedgeMolBonds(c);
        assertEquals( Bond.BondDir.BEGINDASH,m.getBondWithIdx(0).getBondDir() );
        assertEquals( Bond.BondDir.NONE,m.getBondWithIdx(1).getBondDir() );
        assertEquals( Bond.BondDir.NONE,m.getBondWithIdx(2).getBondDir() );
        assertEquals( Bond.BondDir.NONE,m.getBondWithIdx(3).getBondDir() );
    }

    @Test
    public void test2DWithSetAtomLocs() {
        ROMol m = RWMol.MolFromSmiles("C[C@H]1CO1");
        Atom a0 = m.getAtomWithIdx(0);
        Int_Point2D_Map coords = new Int_Point2D_Map();
        coords.set((int)a0.getIdx(), new Point2D(1.0, 1.5));
        RDKFuncs.setPreferCoordGen(false);
        long confIdx = m.compute2DCoords(coords);
        Conformer c = m.getConformer((int) confIdx);
        assertEquals(1.0, c.getAtomPos(a0.getIdx()).getX(), defaultDoubleTol);
        assertEquals(1.5, c.getAtomPos(a0.getIdx()).getY(), defaultDoubleTol);
    }

    @Test
    public void testMatchingDepictions() {
        ROMol template = RWMol.MolFromSmiles("c1nccc2n1ccc2");
        template.compute2DCoords();
        ROMol m = RWMol.MolFromSmiles("c1cccc2ncn3cccc3c21");
        ROMol patt = RWMol.MolFromSmarts("*1****2*1***2");
        Match_Vect mv = m.generateDepictionMatching2DStructure(template,-1,patt);
        assertTrue(mv.size() == 9);
        int[] expected = new int[]{ 6, 5, 4, 12, 11, 7, 8, 9, 10 };
        for (int i = 0; i < mv.size(); ++i) {
            assertTrue(mv.get(i).getFirst() == i);
            assertTrue(mv.get(i).getSecond() == expected[i]);
        }

        // System.out.print(template.MolToMolBlock());
        // System.out.print(m.MolToMolBlock());
        assertEquals(template.MolToMolBlock(), RDKFuncs.MolToMolBlock(template));
        Conformer c1 = template.getConformer();
        Conformer c2 = m.getConformer();
        assertEquals(c1.getAtomPos(0).getX(), c2.getAtomPos(6).getX(), defaultDoubleTol);
        assertEquals(c1.getAtomPos(0).getY(), c2.getAtomPos(6).getY(), defaultDoubleTol);
        assertEquals(c1.getAtomPos(0).getZ(), c2.getAtomPos(6).getZ(), defaultDoubleTol);
    }

    @Test
    public void testGenerateSVG() {
        ROMol m = RWMol.MolFromSmiles("[C@H]1(C)CO1");
        m.compute2DCoords();
        Conformer c = m.getConformer();
        m.WedgeMolBonds(c);
                String svg=m.ToSVG(8,50);
                assertTrue(svg.indexOf("<svg")>-1);
                assertTrue(svg.indexOf("</svg>")>-1);
    }

    @Test
    public void testMolDraw2DSVG() {
        ROMol m = RWMol.MolFromSmiles("[C@H]1(C)CO1");
        m.compute2DCoords();
        Conformer c = m.getConformer();
        m.WedgeMolBonds(c);
        MolDraw2DSVG drawer = new MolDraw2DSVG(300, 300);
        drawer.drawMolecule(m);
        drawer.finishDrawing();
        String svg = drawer.getDrawingText();
        assertTrue(svg.indexOf("<svg") > -1);
        assertTrue(svg.indexOf("</svg>") > -1);
    }

    @Test
    public void testMolDraw2DSVGSingleAtomMol() {
        ROMol m = RWMol.MolFromSmiles("C");
        m.compute2DCoords();
        Conformer c = m.getConformer();
        m.WedgeMolBonds(c);
        MolDraw2DSVG drawer = new MolDraw2DSVG(300, 300);
        drawer.drawMolecule(m);
        drawer.finishDrawing();
        String svg = drawer.getDrawingText();
        assertTrue(svg.indexOf("<svg") > -1);
        assertTrue(svg.indexOf("</svg>") > -1);
    }

    @Test
    public void testNormalizeStraightenDepiction() {
        ROMol noradrenalineMJ = RWMol.MolFromMolBlock("\n" +
"  MJ201100                      \n" +
"\n" +
" 12 12  0  0  1  0  0  0  0  0999 V2000\n" +
"    2.2687    1.0716    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"    1.4437    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"    1.0312    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"    1.4437   -0.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"    0.2062    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"   -0.2062   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"   -1.0312   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"   -1.4437   -1.0716    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"   -1.4437    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"   -2.2687    0.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"   -1.0312    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"   -0.2062    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
"  1  2  1  0  0  0  0\n" +
"  3  2  1  0  0  0  0\n" +
"  3  4  1  6  0  0  0\n" +
"  3  5  1  0  0  0  0\n" +
"  5  6  2  0  0  0  0\n" +
"  6  7  1  0  0  0  0\n" +
"  7  8  1  0  0  0  0\n" +
"  7  9  2  0  0  0  0\n" +
"  9 10  1  0  0  0  0\n" +
"  9 11  1  0  0  0  0\n" +
" 11 12  2  0  0  0  0\n" +
"  5 12  1  0  0  0  0\n" +
"M  END");
        {
            ROMol noradrenalineMJCopy = new RWMol(noradrenalineMJ);
            Conformer conformer0 = noradrenalineMJCopy.getConformer(0);
            Conformer conformer1 = new Conformer(conformer0);
            noradrenalineMJCopy.addConformer(conformer1, true);
            assertTrue(noradrenalineMJ.calcRMS(noradrenalineMJCopy, 0, 0) < 1.e-5);
            assertTrue(noradrenalineMJ.calcRMS(noradrenalineMJCopy, 0, 1) < 1.e-5);
            double scalingFactor = noradrenalineMJCopy.normalizeDepiction(1);
            assertTrue(noradrenalineMJ.calcRMS(noradrenalineMJCopy, 0, 0) < 1.e-5);
            assertTrue(noradrenalineMJ.calcRMS(noradrenalineMJCopy, 0, 1) > 1.e-5);
            assertEquals(scalingFactor, 1.875, 1.e-3);
            Conformer conformer2 = new Conformer(conformer1);
            noradrenalineMJCopy.addConformer(conformer2, true);
            Point3D bond10_11Conf0 = conformer0.getAtomPos(11).minus(conformer0.getAtomPos(10));
            assertEquals(bond10_11Conf0.getX(), 0.825, 1.e-3);
            assertEquals(bond10_11Conf0.getY(), 0.0, 1.e-3);
            Point3D bond10_11Conf1 = conformer1.getAtomPos(11).minus(conformer1.getAtomPos(10));
            assertEquals(bond10_11Conf1.getX(), 1.513, 1.e-3);
            assertEquals(bond10_11Conf1.getY(), -0.321, 1.e-3);
            noradrenalineMJCopy.straightenDepiction(1);
            bond10_11Conf1 = conformer1.getAtomPos(11).minus(conformer1.getAtomPos(10));
            assertEquals(bond10_11Conf1.getX(), 1.340, 1.e-3);
            assertEquals(bond10_11Conf1.getY(), -0.773, 1.e-3);
            Point3D bond4_11Conf1 = conformer1.getAtomPos(11).minus(conformer1.getAtomPos(4));
            assertEquals(bond4_11Conf1.getX(), 0.0, 1.e-3);
            assertEquals(bond4_11Conf1.getY(), 1.547, 1.e-3);
            noradrenalineMJCopy.straightenDepiction(2, true);
            Point3D bond10_11Conf2 = conformer2.getAtomPos(11).minus(conformer2.getAtomPos(10));
            assertEquals(bond10_11Conf2.getX(), 1.547, 1.e-3);
            assertEquals(bond10_11Conf2.getY(), 0.0, 1.e-3);
            Point3D bond4_11Conf2 = conformer2.getAtomPos(11).minus(conformer2.getAtomPos(4));
            assertEquals(bond4_11Conf2.getX(), -0.773, 1.e-3);
            assertEquals(bond4_11Conf2.getY(), 1.339, 1.e-3);
            noradrenalineMJCopy.delete();
            conformer1.delete();
            conformer2.delete();
        }
        {
            ROMol noradrenalineMJCopy = new RWMol(noradrenalineMJ);
            Conformer conformer0 = noradrenalineMJCopy.getConformer(0);
            Conformer conformer1 = new Conformer(conformer0);
            noradrenalineMJCopy.addConformer(conformer1, true);
            double scalingFactor = noradrenalineMJCopy.normalizeDepiction(1, -1);
            assertTrue(noradrenalineMJ.calcRMS(noradrenalineMJCopy, 0, 0) < 1.e-5);
            assertTrue(noradrenalineMJ.calcRMS(noradrenalineMJCopy, 0, 1) > 1.e-5);
            assertEquals(scalingFactor, 1.875, 1.e-3);
            Conformer conformer2 = new Conformer(conformer1);
            noradrenalineMJCopy.addConformer(conformer2, true);
            Point3D bond10_11Conf0 = conformer0.getAtomPos(11).minus(conformer0.getAtomPos(10));
            assertEquals(bond10_11Conf0.getX(), 0.825, 1.e-3);
            assertEquals(bond10_11Conf0.getY(), 0.0, 1.e-3);
            Point3D bond10_11Conf1 = conformer1.getAtomPos(11).minus(conformer1.getAtomPos(10));
            assertEquals(bond10_11Conf1.getX(), 0.321, 1.e-3);
            assertEquals(bond10_11Conf1.getY(), 1.513, 1.e-3);
            noradrenalineMJCopy.straightenDepiction(1);
            bond10_11Conf1 = conformer1.getAtomPos(11).minus(conformer1.getAtomPos(10));
            assertEquals(bond10_11Conf1.getX(), 0.0, 1.e-3);
            assertEquals(bond10_11Conf1.getY(), 1.547, 1.e-3);
            noradrenalineMJCopy.straightenDepiction(2, true);
            Point3D bond10_11Conf2 = conformer2.getAtomPos(11).minus(conformer2.getAtomPos(10));
            assertEquals(bond10_11Conf2.getX(), bond10_11Conf1.getX(), 1.e-3);
            assertEquals(bond10_11Conf2.getY(), bond10_11Conf1.getY(), 1.e-3);
            noradrenalineMJCopy.delete();
            conformer1.delete();
            conformer2.delete();
        }
        {
            ROMol noradrenalineMJCopy = new RWMol(noradrenalineMJ);
            Conformer conformer0 = noradrenalineMJCopy.getConformer(0);
            Conformer conformer1 = new Conformer(conformer0);
            noradrenalineMJCopy.addConformer(conformer1, true);
            double scalingFactor = noradrenalineMJCopy.normalizeDepiction(1, 0, 3.0);
            assertTrue(noradrenalineMJ.calcRMS(noradrenalineMJCopy, 0, 0) < 1.e-5);
            assertTrue(noradrenalineMJ.calcRMS(noradrenalineMJCopy, 0, 1) > 1.e-5);
            assertEquals(scalingFactor, 3.0, 1.e-3);
            Conformer conformer2 = new Conformer(conformer1);
            noradrenalineMJCopy.addConformer(conformer2, true);
            Conformer conformer3 = new Conformer(conformer1);
            noradrenalineMJCopy.addConformer(conformer3, true);
            Point3D bond10_11Conf0 = conformer0.getAtomPos(11).minus(conformer0.getAtomPos(10));
            assertEquals(bond10_11Conf0.getX(), 0.825, 1.e-3);
            assertEquals(bond10_11Conf0.getY(), 0.0, 1.e-3);
            Point3D bond10_11Conf1 = conformer1.getAtomPos(11).minus(conformer1.getAtomPos(10));
            assertEquals(bond10_11Conf1.getX(), 2.475, 1.e-3);
            assertEquals(bond10_11Conf1.getY(), 0.0, 1.e-3);
            noradrenalineMJCopy.straightenDepiction(1);
            bond10_11Conf1 = conformer1.getAtomPos(11).minus(conformer1.getAtomPos(10));
            assertEquals(bond10_11Conf1.getX(), 2.143, 1.e-3);
            assertEquals(bond10_11Conf1.getY(), -1.237, 1.e-3);
            Point3D bond4_11Conf1 = conformer1.getAtomPos(11).minus(conformer1.getAtomPos(4));
            assertEquals(bond4_11Conf1.getX(), 0.0, 1.e-3);
            assertEquals(bond4_11Conf1.getY(), 2.475, 1.e-3);
            noradrenalineMJCopy.straightenDepiction(2, true);
            Point3D bond10_11Conf2 = conformer2.getAtomPos(11).minus(conformer2.getAtomPos(10));
            Point3D bond10_11Conf3 = conformer3.getAtomPos(11).minus(conformer3.getAtomPos(10));
            assertEquals(bond10_11Conf2.getX(), bond10_11Conf3.getX(), 1.e-3);
            assertEquals(bond10_11Conf2.getY(), bond10_11Conf3.getY(), 1.e-3);
            Point3D bond4_11Conf2 = conformer2.getAtomPos(11).minus(conformer2.getAtomPos(4));
            Point3D bond4_11Conf3 = conformer3.getAtomPos(11).minus(conformer3.getAtomPos(4));
            assertEquals(bond4_11Conf2.getX(), bond4_11Conf3.getX(), 1.e-3);
            assertEquals(bond4_11Conf2.getY(), bond4_11Conf3.getY(), 1.e-3);
            noradrenalineMJCopy.delete();
            conformer1.delete();
            conformer2.delete();
            conformer3.delete();
        }
    }

    @Test
    public void testGetBestRMS() {
        String rdpath = System.getenv("RDBASE");
        if (rdpath == null)
            org.junit.Assert.fail("No definition for RDBASE");
        File base = new File(rdpath);
        File testFile = new File(base, "Code" + File.separator + "GraphMol"
                + File.separator + "MolAlign" + File.separator + "test_data"
                + File.separator + "probe_mol.sdf");
        String fName = testFile.getAbsolutePath();
        SDMolSupplier supplier = new SDMolSupplier(fName, true, false);
        supplier.next();
        ROMol prb = supplier.next();
        ROMol ref = supplier.next();
        ROMol prbCopy1 = new ROMol(prb);
        ROMol prbCopy2 = new ROMol(prb);
        ROMol prbCopy3 = new ROMol(prb);
        Transform3D bestTrans = new Transform3D();
        Match_Vect bestMatch = new Match_Vect();

        // alignMol() would return this for the rms: 2.50561
        // But the best rms is: 2.43449
        double rmsdInPlace = prbCopy1.calcRMS(ref);
        assertEquals(rmsdInPlace, 2.6026, 1.e-3);
        double rmsd = prb.getBestRMS(ref);
        assertEquals(rmsd, 2.43449, 1.e-3);
        double rmsdCopy = prbCopy1.getBestAlignmentTransform(ref, bestTrans, bestMatch);
        assertEquals(rmsd, rmsdCopy, 1.e-3);
        assertEquals(bestMatch.size(), ref.getNumAtoms());

        SmilesParserParams params = new SmilesParserParams();
        params.setRemoveHs(false);
        RWMol scaffold = RWMol.MolFromSmiles("N1C([H])([H])C([H])([H])C([H])([H])[N+]([H])([H])C([H])([H])C1([H])[H]", params);
        Match_Vect scaffoldMatch = ref.getSubstructMatch(scaffold);
        assertFalse(scaffoldMatch.isEmpty());
        int nAtoms = (int)ref.getNumAtoms();
        boolean[] scaffoldIndices = new boolean[nAtoms];
        Arrays.fill(scaffoldIndices, false);
        for (int i = 0; i < scaffoldMatch.size(); ++i) {
            int pairSecond = scaffoldMatch.get(i).getSecond();
            scaffoldIndices[pairSecond] = true;
        }
        Match_Vect_Vect matches = ref.getSubstructMatches(prb, false);
        Match_Vect_Vect matchesPruned = new Match_Vect_Vect(matches.size());
        for (int i = 0; i < matches.size(); ++i) {
            Match_Vect match = matches.get(i);
            Match_Vect matchPruned = new Match_Vect();
            for (int j = 0; j < match.size(); ++j) {
                int pairSecond = match.get(j).getSecond();
                if (scaffoldIndices[pairSecond]) {
                    matchPruned.add(match.get(j));
                }
            }
            matchesPruned.set(i, matchPruned);
            matchPruned.delete();
        }
        rmsdInPlace = prbCopy2.calcRMS(ref, -1, -1, matchesPruned);
        assertEquals(rmsdInPlace, 2.5672, 1.e-3);
        rmsd = prb.getBestRMS(ref, -1, -1, matchesPruned);
        assertEquals(rmsd, 1.14329, 1.e-3);
        rmsdCopy = prbCopy2.getBestAlignmentTransform(ref, bestTrans, bestMatch, -1, -1, matchesPruned);
        assertEquals(rmsd, rmsdCopy, 1.e-3);
        assertEquals(bestMatch.size(), scaffoldMatch.size());
        DoubleVector weights = new DoubleVector(nAtoms, 1.0);
        for (int i = 0; i < nAtoms; ++i) {
            if (scaffoldIndices[i]) {
                weights.setVal(i, 100.0);
            }
        }
        rmsdInPlace = prbCopy3.calcRMS(ref, -1, -1, matches, 1000, true, weights);
        assertEquals(rmsdInPlace, 17.7959, 1.e-3);
        rmsd = prb.getBestRMS(ref, -1, -1, matches, 1000, true, weights);
        assertEquals(rmsd, 10.9681, 1.e-3);
        rmsdCopy = prbCopy3.getBestAlignmentTransform(ref, bestTrans, bestMatch, -1, -1, matches, 1000, true, weights);
        assertEquals(rmsd, rmsdCopy, 1.e-3);
        assertEquals(bestMatch.size(), ref.getNumAtoms());
        supplier.delete();
        prbCopy1.delete();
        prbCopy2.delete();
        prbCopy3.delete();
        params.delete();
        matchesPruned.delete();
        weights.delete();
    }

    @Test
    public void testPrepareMolForDrawing() {
        RWMol m = RWMol.MolFromSmiles("c1ccccc1");
        RDKFuncs.prepareMolForDrawing(m);
        assertEquals(m.getNumConformers(), 1);
        assertTrue(m.getBondBetweenAtoms(0, 1).getBondType() !=
                    Bond.BondType.AROMATIC);
        assertTrue(m.getBondBetweenAtoms(0, 1).getIsAromatic());

        MolDraw2DSVG drawer = new MolDraw2DSVG(300, 300);
        drawer.drawMolecule(m);
        drawer.finishDrawing();
        String svg = drawer.getDrawingText();
        assertTrue(svg.indexOf("<svg") > -1);
        assertTrue(svg.indexOf("</svg>") > -1);
    }

    @Test
    public void testMolDraw2DHighlight() {
        RWMol m = RWMol.MolFromSmiles("CCCCCOC");
        RDKFuncs.prepareMolForDrawing(m);
        Int_Vect hats = new Int_Vect();
        hats.add(0);
        hats.add(1);
        hats.add(2);

        Int_Vect hbs = new Int_Vect();
        hbs.add(0);
        hbs.add(1);
        hbs.add(2);

        ColourPalette atCs = new ColourPalette();
        atCs.set(0, new DrawColour(1, 1, 0));
        atCs.set(1, new DrawColour(1, 0, 1));
        atCs.set(2, new DrawColour(0, 1, 1));
        ColourPalette bCs = new ColourPalette();

        MolDraw2DSVG drawer = new MolDraw2DSVG(300, 300);
        drawer.drawMolecule(m, "THE_LEGEND", hats, hbs, atCs, bCs);
        drawer.finishDrawing();
        String svg = drawer.getDrawingText();
        // System.out.print(svg);
        assertTrue(svg.indexOf("<svg") > -1);
        assertTrue(svg.indexOf("</svg>") > -1);
        assertTrue(svg.indexOf("fill:#FFFF00;") > -1);
        assertTrue(svg.indexOf("fill:#FF00FF;") > -1);
        assertTrue(svg.indexOf("fill:#00FFFF;") > -1);
        // default line color:
        assertTrue(svg.indexOf("stroke:#FF7F7F;") == -1);
        assertTrue(svg.indexOf("stroke:#FFFF00;") > -1);
        assertTrue(svg.indexOf("stroke:#FF00FF;") > -1);
    }

    @Test
    public void testMolDraw2DContours() {
        // really just a test to make sure this can be called
        RWMol m = RWMol.MolFromSmiles("CCCCCOC");
        RDKFuncs.prepareMolForDrawing(m);
        Point2D_Vect cents = new Point2D_Vect();
        Double_Vect weights = new Double_Vect();
        Double_Vect widths = new Double_Vect();
        for(int i=0;i<m.getNumAtoms();++i){
            cents.add(new Point2D(m.getConformer().getAtomPos(i)));
            weights.add(1.0);
            widths.add(0.4);
        }
        MolDraw2DSVG drawer = new MolDraw2DSVG(300, 300);
        drawer.clearDrawing();
        RDKFuncs.ContourAndDrawGaussians(drawer,cents,weights,widths,10);
        drawer.drawOptions().setClearBackground(false);
        drawer.drawMolecule(m);
        drawer.finishDrawing();
        String svg = drawer.getDrawingText();
        System.out.print(svg);
    }

    @Test
    public void testStrictParsingAndLogging() {
        RDKFuncs.InitLogs();
        String badMolBlock = "\n" +
            "  MJ201100                      \n" +
            "\n" +
            "  3  2  0  0  0  0  0  0  0  0999 V2000\n" +
            "   -0.3572   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    0.3572    0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "    1.0717    0.6187    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
            "  1  2  1  0  0  0  0\n" +
            "  2  3  3  0  0  0  0\n" +
            "M  STY  1   1 SUP\n" +
            "M  SAL   1  2   2   3\n" +
            "M  SMT   1 CN\n" +
            "M  SBL   1  1   1\n" +
            "M  SAP   1  1   2\n" +
            "M  END\n";
        boolean exceptionThrown = false;
        boolean molIsValid = false;
        ROMol mol = null;
        try {
            mol = RWMol.MolFromMolBlock(badMolBlock);
        } catch(Exception e) {
            exceptionThrown = true;
        } finally {
            if (mol != null) {
                molIsValid = true;
                mol.delete();
            }
        }
        assertTrue(exceptionThrown);
        assertFalse(molIsValid);
        exceptionThrown = false;
        BufferedReader reader = null;
        try {
            String filename = "java_warning_log.txt";
            RDKFuncs.getRdWarningLog().SetTee(filename);
            mol = RWMol.MolFromMolBlock(badMolBlock, true, true, false);
            reader = new BufferedReader(new FileReader(filename));
            String line = reader.readLine();
            assertTrue(line != null);
            assertTrue(line.contains("SGroup SAP line too short"));
        } catch(Exception e) {
            exceptionThrown = true;
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e) {}
            }
            if (mol != null) {
                molIsValid = true;
                mol.delete();
            }
        }
        assertFalse(exceptionThrown);
        assertTrue(molIsValid);
    }

    @Test
    public void testMostSubstitutedCoreMatch() {
        RWMol core = RWMol.MolFromSmarts("[*:1]c1cc([*:2])ccc1[*:3]");
        RWMol orthoMeta = RWMol.MolFromSmiles("c1ccc(-c2ccc(-c3ccccc3)c(-c3ccccc3)c2)cc1");
        RWMol ortho = RWMol.MolFromSmiles("c1ccc(-c2ccccc2-c2ccccc2)cc1");
        RWMol meta = RWMol.MolFromSmiles("c1ccc(-c2cccc(-c3ccccc3)c2)cc1");
        RWMol biphenyl = RWMol.MolFromSmiles("c1ccccc1-c1ccccc1");
        RWMol phenyl = RWMol.MolFromSmiles("c1ccccc1");

        class NumHsMatchingDummies {
            public int get(RWMol mol, RWMol core, Match_Vect match) {
                int count = 0;
                for (int i = 0; i < match.size(); ++i) {
                    Int_Pair pair = match.get(i);
                    if (core.getAtomWithIdx(pair.getFirst()).getAtomicNum() == 0 &&
                        mol.getAtomWithIdx(pair.getSecond()).getAtomicNum() == 1) {
                        ++count;
                    }
                }
                return count;
            }
        };

        NumHsMatchingDummies numHsMatchingDummies = new NumHsMatchingDummies();
        RWMol_Vect mols = new RWMol_Vect(5);
        mols.set(0, orthoMeta);
        mols.set(1, ortho);
        mols.set(2, meta);
        mols.set(3, biphenyl);
        mols.set(4, phenyl);
        int[] expected = new int[]{ 0, 1, 1, 2, 3 };
        assertTrue(mols.size() == expected.length);
        for (int i = 0; i < expected.length; ++i) {
            RWMol mol = mols.get(i);
            int res = expected[i];
            RDKFuncs.addHs(mol);
            Match_Vect_Vect matches = mol.getSubstructMatches(core);
            Match_Vect bestMatch = mol.getMostSubstitutedCoreMatch(core, matches);
            assertTrue(numHsMatchingDummies.get(mol, core, bestMatch) == res);
            int[] ctrlCounts = new int[(int)matches.size()];
            for (int j = 0; j < ctrlCounts.length; ++j) {
                ctrlCounts[j] = numHsMatchingDummies.get(mol, core, matches.get(j));
            }
            Arrays.sort(ctrlCounts);
            int[] sortedCounts = new int[ctrlCounts.length];
            Match_Vect_Vect sortedMatches = mol.sortMatchesByDegreeOfCoreSubstitution(core, matches);
            for (int j = 0; j < sortedMatches.size(); ++j) {
                sortedCounts[j] = numHsMatchingDummies.get(mol, core, sortedMatches.get(j));
            }
            assertTrue(Arrays.equals(ctrlCounts, sortedCounts));
            matches.delete();
            sortedMatches.delete();
        }
        Match_Vect_Vect emptyMatches = new Match_Vect_Vect();
        boolean raised = false;
        try {
            orthoMeta.getMostSubstitutedCoreMatch(core, emptyMatches);
        } catch (Exception e) {
            raised = true;
        }
        assertTrue(raised);
        raised = false;
        try {
            orthoMeta.sortMatchesByDegreeOfCoreSubstitution(core, emptyMatches);
        } catch (Exception e) {
            raised = true;
        }
        assertTrue(raised);
        emptyMatches.delete();
        mols.delete();
        core.delete();
        orthoMeta.delete();
        ortho.delete();
        meta.delete();
        biphenyl.delete();
        phenyl.delete();
    }

    @Test
    public void testStereoChemFunctions() {
        boolean useLegacyStereo = RDKFuncs.getUseLegacyStereoPerception();
        boolean allowNonTetrahedralChirality = RDKFuncs.getAllowNontetrahedralChirality();
        try {
            RWMol m = RWMol.MolFromSmiles("CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wU:1.0|");
            assertTrue(m != null);
            assertEquals(m.getBondWithIdx(0).getProp("_MolFileBondCfg"), "1");
            assertEquals(m.getAtomWithIdx(1).getChiralTag(), Atom.ChiralType.CHI_TETRAHEDRAL_CW);
            m.delete();
            m = RWMol.MolFromSmiles("CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wD:1.0|");
            assertTrue(m != null);
            assertEquals(m.getBondWithIdx(0).getProp("_MolFileBondCfg"), "3");
            assertEquals(m.getAtomWithIdx(1).getChiralTag(), Atom.ChiralType.CHI_TETRAHEDRAL_CCW);
            m.invertMolBlockWedgingInfo();
            assertEquals(m.getBondWithIdx(0).getProp("_MolFileBondCfg"), "1");
            m.reapplyMolBlockWedging();
            RDKFuncs.assignChiralTypesFromBondDirs(m);
            assertEquals(m.getAtomWithIdx(1).getChiralTag(), Atom.ChiralType.CHI_TETRAHEDRAL_CW);
            m.invertMolBlockWedgingInfo();
            assertEquals(m.getBondWithIdx(0).getProp("_MolFileBondCfg"), "3");
            m.reapplyMolBlockWedging();
            RDKFuncs.assignChiralTypesFromBondDirs(m);
            assertEquals(m.getAtomWithIdx(1).getChiralTag(), Atom.ChiralType.CHI_TETRAHEDRAL_CCW);
            m.clearMolBlockWedgingInfo();
            m.getAtomWithIdx(1).setChiralTag(Atom.ChiralType.CHI_UNSPECIFIED);
            assertFalse(m.getBondWithIdx(0).hasProp("_MolFileBondCfg"));
            m.reapplyMolBlockWedging();
            RDKFuncs.assignChiralTypesFromBondDirs(m);
            assertEquals(m.getAtomWithIdx(1).getChiralTag(), Atom.ChiralType.CHI_UNSPECIFIED);
            m.delete();
            RDKFuncs.setUseLegacyStereoPerception(true);
            m = RWMol.MolFromSmiles("O[C@@]1(C)C/C(/C1)=C(/C)\\CC");
            assertTrue(m != null);
            assertEquals(m.MolToSmiles(), "CCC(C)=C1CC(C)(O)C1");
            m.delete();
            RDKFuncs.setUseLegacyStereoPerception(false);
            m = RWMol.MolFromSmiles("O[C@@]1(C)C/C(/C1)=C(/C)\\CC");
            assertTrue(m != null);
            assertEquals(m.MolToSmiles(), "CC/C(C)=C1\\C[C@](C)(O)C1");
            m.delete();
            RDKFuncs.setUseLegacyStereoPerception(useLegacyStereo);
            String ctab = "\n" +
                "  Mrv2108 09132105183D          \n" +
                "\n" +
                "  5  4  0  0  0  0            999 V2000\n" +
                "   -1.2500    1.4518    0.0000 Pt  0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -1.2500    2.2768    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -0.4250    1.4518    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -2.0750    1.4518    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -1.2500    0.6268    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  1  0  0  0  0\n" +
                "  1  3  1  0  0  0  0\n" +
                "  1  4  1  0  0  0  0\n" +
                "  1  5  1  0  0  0  0\n" +
                "M  END\n";
            RDKFuncs.setAllowNontetrahedralChirality(true);
            m = RWMol.MolFromMolBlock(ctab);
            assertTrue(m != null);
            assertEquals(m.MolToSmiles(), "F[Pt@SP3](F)(Cl)Cl");
            m.delete();
            RDKFuncs.setAllowNontetrahedralChirality(false);
            m = RWMol.MolFromMolBlock(ctab);
            assertTrue(m != null);
            assertEquals(m.MolToSmiles(), "F[Pt](F)(Cl)Cl");
            m.delete();
            RDKFuncs.setAllowNontetrahedralChirality(allowNonTetrahedralChirality);
        } finally {
            RDKFuncs.setUseLegacyStereoPerception(useLegacyStereo);
            RDKFuncs.setAllowNontetrahedralChirality(allowNonTetrahedralChirality);
        }
    }
    
    public static void main(String args[]) {
        org.junit.runner.JUnitCore.main("org.RDKit.Chemv2Tests");
    }
}
