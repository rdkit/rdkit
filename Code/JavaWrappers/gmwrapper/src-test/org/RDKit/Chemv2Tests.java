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
	public void testEditingPersisting()	{
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
		coords.set((int) a0.getIdx(), new Point2D(1.0, 1.5));
		RDKFuncs.setPreferCoordGen(false);
		long confIdx = m.compute2DCoords(coords);
		Conformer c = m.getConformer((int) confIdx);
		assertEquals(1.0, c.getAtomPos(a0.getIdx()).getX(), defaultDoubleTol);
		assertEquals(1.5, c.getAtomPos(a0.getIdx()).getY(), defaultDoubleTol);
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
          assertTrue(svg.indexOf("THE_LEGEND") > -1);
          assertTrue(svg.indexOf("fill:#FFFF00;") > -1);
          assertTrue(svg.indexOf("fill:#FF00FF;") > -1);
          assertTrue(svg.indexOf("fill:#00FFFF;") > -1);
          // default line color:
          assertTrue(svg.indexOf("stroke:#FF7F7F;") > -1);
        }

        public static void main(String args[]) {
		org.junit.runner.JUnitCore.main("org.RDKit.Chemv2Tests");
	}

}
